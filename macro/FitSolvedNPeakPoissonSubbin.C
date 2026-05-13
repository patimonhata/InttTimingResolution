#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "Utility/sPhenixStyle.C"

void PutCaption(double first_margin = -0.1, double pos_x = 0.2, double pos_y = 0.9);

namespace {

struct MultiplicityConfig {
  const char* dir_suffix;
  const char* label;
  int color;
};

struct FitConfig {
  std::string group_label = "Group2";
  std::string base_output_dir =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
  std::string output_dir =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
  std::string histogram_name = "h_solved_n";
  int origin_bin = 8;
  int fit_first_bin = 9;
  int fit_last_bin = 13;
  int subbin_factor = 10;
  bool draw_full_x_range = true;
};

struct FitSummary {
  std::string label;
  int status = -1;
  double amplitude = std::numeric_limits<double>::quiet_NaN();
  double amplitude_error = std::numeric_limits<double>::quiet_NaN();
  double lambda = std::numeric_limits<double>::quiet_NaN();
  double lambda_error = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = 0;
};

double gModelOriginX = 0.0;
double gModelCoarseBinWidth = 1.0;
int gModelSubbinFactor = 10;

std::vector<MultiplicityConfig> BuildMultiplicityConfigs() {
  return {
      {"mul10to100", "Low multiplicity", kRed + 1},
      {"mul100to200", "Middle multiplicity", kBlue + 1},
      {"mul200to1000", "High multiplicity", kGreen + 2},
  };
}

double SumPoissonTerms(int first_k, int last_k, double lambda) {
  if (lambda <= 0.0 || last_k < first_k) {
    return 0.0;
  }

  double sum = 0.0;
  for (int k = first_k; k <= last_k; ++k) {
    sum += TMath::Poisson(k, lambda);
  }
  return sum;
}

double CoarseBinExpectation(double* x, double* par) {
  const double amplitude = par[0];
  const double lambda = par[1];
  if (amplitude < 0.0 || lambda <= 0.0 || gModelSubbinFactor <= 0 ||
      gModelCoarseBinWidth <= 0.0) {
    return 0.0;
  }

  const double coarse_position = (x[0] - gModelOriginX) / gModelCoarseBinWidth;
  if (coarse_position < 0.0) {
    return 0.0;
  }

  const int coarse_index = static_cast<int>(std::floor(coarse_position));
  const int first_k = coarse_index * gModelSubbinFactor;
  const int last_k = first_k + gModelSubbinFactor - 1;
  return amplitude * SumPoissonTerms(first_k, last_k, lambda);
}

TH1D* LoadSolvedHistogram(const FitConfig& config, const MultiplicityConfig& multiplicity) {
  const TString input_path = Form("%s/shifted_sum_reconstruction_%s_%s/solved_n.root",
                                  config.base_output_dir.c_str(),
                                  config.group_label.c_str(),
                                  multiplicity.dir_suffix);

  TFile* input_file = TFile::Open(input_path, "READ");
  if (input_file == nullptr || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << input_path << std::endl;
    delete input_file;
    return nullptr;
  }

  TH1D* input_hist =
      dynamic_cast<TH1D*>(input_file->Get(config.histogram_name.c_str()));
  if (input_hist == nullptr) {
    std::cerr << "ERROR: " << config.histogram_name << " is missing in "
              << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(
      input_hist->Clone(Form("%s_%s_%s",
                             config.histogram_name.c_str(),
                             config.group_label.c_str(),
                             multiplicity.dir_suffix)));
  if (hist == nullptr) {
    std::cerr << "ERROR: Failed to clone " << config.histogram_name << " from "
              << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  hist->SetDirectory(nullptr);
  input_file->Close();
  delete input_file;
  return hist;
}

bool ValidateConfig(const FitConfig& config, const TH1D& hist) {
  if (config.subbin_factor <= 0) {
    std::cerr << "ERROR: subbin_factor must be positive." << std::endl;
    return false;
  }

  if (config.origin_bin < 1 || config.origin_bin > hist.GetNbinsX()) {
    std::cerr << "ERROR: origin_bin = " << config.origin_bin
              << " is outside histogram range." << std::endl;
    return false;
  }

  if (config.fit_first_bin < 1 || config.fit_last_bin > hist.GetNbinsX() ||
      config.fit_first_bin > config.fit_last_bin) {
    std::cerr << "ERROR: Invalid fit bin range [" << config.fit_first_bin << ", "
              << config.fit_last_bin << "]." << std::endl;
    return false;
  }

  if (config.fit_first_bin <= config.origin_bin) {
    std::cerr << "ERROR: fit_first_bin must be larger than origin_bin."
              << std::endl;
    return false;
  }

  const double first_width = hist.GetXaxis()->GetBinWidth(1);
  for (int bin = 2; bin <= hist.GetNbinsX(); ++bin) {
    if (std::abs(hist.GetXaxis()->GetBinWidth(bin) - first_width) > 1e-9) {
      std::cerr << "ERROR: Non-uniform bin width is not supported." << std::endl;
      return false;
    }
  }

  return true;
}

double EstimateInitialLambda(const TH1D& hist, const FitConfig& config) {
  int peak_bin = config.fit_first_bin;
  double peak_value = hist.GetBinContent(config.fit_first_bin);
  for (int bin = config.fit_first_bin + 1; bin <= config.fit_last_bin; ++bin) {
    if (hist.GetBinContent(bin) > peak_value) {
      peak_value = hist.GetBinContent(bin);
      peak_bin = bin;
    }
  }

  const int coarse_index = peak_bin - config.origin_bin;
  const double coarse_center =
      (static_cast<double>(coarse_index) + 0.5) * config.subbin_factor;
  return std::max(coarse_center - 0.5, 1.0);
}

double EstimateInitialAmplitude(const TH1D& hist, const FitConfig& config, double lambda) {
  double observed_sum = 0.0;
  for (int bin = config.fit_first_bin; bin <= config.fit_last_bin; ++bin) {
    observed_sum += std::max(0.0, hist.GetBinContent(bin));
  }

  double probability_sum = 0.0;
  for (int bin = config.fit_first_bin; bin <= config.fit_last_bin; ++bin) {
    const int coarse_index = bin - config.origin_bin;
    const int first_k = coarse_index * config.subbin_factor;
    const int last_k = first_k + config.subbin_factor - 1;
    probability_sum += SumPoissonTerms(first_k, last_k, lambda);
  }

  if (probability_sum <= 0.0) {
    return std::max(observed_sum, 1.0);
  }
  return std::max(observed_sum / probability_sum, 1.0);
}

TF1* BuildFitFunction(const TH1D& hist,
                      const FitConfig& config,
                      const MultiplicityConfig& multiplicity,
                      const char* suffix) {
  const double x_min = hist.GetXaxis()->GetXmin();
  const double x_max = hist.GetXaxis()->GetXmax();
  const TString function_name =
      Form("f_poisson_subbin_%s_%s_%s",
           config.group_label.c_str(),
           multiplicity.dir_suffix,
           suffix);

  return new TF1(function_name, CoarseBinExpectation, x_min, x_max, 2);
}

TH1D* BuildFineModelHistogram(const TH1D& coarse_hist,
                              const FitConfig& config,
                              const MultiplicityConfig& multiplicity,
                              double amplitude,
                              double lambda) {
  const int fine_bins = coarse_hist.GetNbinsX() * config.subbin_factor;
  const double x_min = coarse_hist.GetXaxis()->GetXmin();
  const double x_max = coarse_hist.GetXaxis()->GetXmax();

  TH1D* fine_hist = new TH1D(Form("h_fine_model_%s_%s",
                                  config.group_label.c_str(),
                                  multiplicity.dir_suffix),
                             "",
                             fine_bins,
                             x_min,
                             x_max);
  fine_hist->SetDirectory(nullptr);

  for (int fine_bin = 1; fine_bin <= fine_bins; ++fine_bin) {
    const double x = fine_hist->GetBinCenter(fine_bin);
    const double fine_position =
        (x - gModelOriginX) / (gModelCoarseBinWidth / config.subbin_factor);
    if (fine_position < 0.0) {
      fine_hist->SetBinContent(fine_bin, 0.0);
      continue;
    }

    const int k = static_cast<int>(std::floor(fine_position));
    fine_hist->SetBinContent(fine_bin, amplitude * TMath::Poisson(k, lambda));
  }

  fine_hist->SetLineColor(multiplicity.color);
  fine_hist->SetLineStyle(7);
  fine_hist->SetLineWidth(2);
  return fine_hist;
}

FitSummary FitSingleHistogram(TH1D& hist,
                              const FitConfig& config,
                              const MultiplicityConfig& multiplicity,
                              TF1*& draw_function,
                              TH1D*& fine_model_hist) {
  FitSummary summary;
  summary.label = multiplicity.label;

  gModelOriginX = hist.GetXaxis()->GetBinLowEdge(config.origin_bin);
  gModelCoarseBinWidth = hist.GetXaxis()->GetBinWidth(config.origin_bin);
  gModelSubbinFactor = config.subbin_factor;

  const double initial_lambda = EstimateInitialLambda(hist, config);
  const double initial_amplitude = EstimateInitialAmplitude(hist, config, initial_lambda);

  std::unique_ptr<TF1> fit_function(BuildFitFunction(hist, config, multiplicity, "fit"));
  fit_function->SetRange(hist.GetXaxis()->GetBinLowEdge(config.fit_first_bin),
                         hist.GetXaxis()->GetBinUpEdge(config.fit_last_bin));
  fit_function->SetParName(0, "A");
  fit_function->SetParName(1, "lambda");
  fit_function->SetParameter(0, initial_amplitude);
  fit_function->SetParameter(1, initial_lambda);
  fit_function->SetParLimits(0, 0.0, 1.0e9);
  fit_function->SetParLimits(1, 0.01, 1.0e3);
  fit_function->SetLineColor(multiplicity.color);
  fit_function->SetLineWidth(2);

  TFitResultPtr fit_result = hist.Fit(fit_function.get(), "RS0QWW");
  summary.status = fit_result.Get() ? fit_result->Status() : -1;
  if (fit_result.Get() == nullptr || summary.status != 0) {
    std::cerr << "ERROR: Fit failed for " << multiplicity.label
              << " with status " << summary.status << std::endl;
    return summary;
  }

  summary.amplitude = fit_function->GetParameter(0);
  summary.amplitude_error = fit_function->GetParError(0);
  summary.lambda = fit_function->GetParameter(1);
  summary.lambda_error = fit_function->GetParError(1);
  summary.chi2 = fit_function->GetChisquare();
  summary.ndf = fit_function->GetNDF();

  draw_function = BuildFitFunction(hist, config, multiplicity, "draw");
  draw_function->SetParameters(summary.amplitude, summary.lambda);
  draw_function->SetLineColor(multiplicity.color);
  draw_function->SetLineWidth(2);

  fine_model_hist = BuildFineModelHistogram(
      hist, config, multiplicity, summary.amplitude, summary.lambda);
  return summary;
}

void DrawFitGuides(const TH1D& hist, const FitConfig& config, double y_max) {
  const double origin_x = hist.GetXaxis()->GetBinLowEdge(config.origin_bin);
  const double fit_x_min = hist.GetXaxis()->GetBinLowEdge(config.fit_first_bin);
  const double fit_x_max = hist.GetXaxis()->GetBinUpEdge(config.fit_last_bin);

  TLine origin_line(origin_x, 0.0, origin_x, y_max);
  origin_line.SetLineStyle(2);
  origin_line.SetLineColor(kGray + 2);
  origin_line.Draw();

  TLine left_line(fit_x_min, 0.0, fit_x_min, y_max);
  TLine right_line(fit_x_max, 0.0, fit_x_max, y_max);
  left_line.SetLineStyle(3);
  right_line.SetLineStyle(3);
  left_line.SetLineColor(kBlack);
  right_line.SetLineColor(kBlack);
  left_line.Draw();
  right_line.Draw();
}

void DrawSinglePad(TH1D& hist,
                   const FitConfig& config,
                   const MultiplicityConfig& multiplicity,
                   const FitSummary& summary,
                   TF1* draw_function,
                   TH1D* fine_model_hist,
                   bool draw_legend) {
  hist.SetStats(0);
  hist.SetTitle("");
  hist.SetLineColor(kBlack);
  hist.SetLineWidth(2);
  hist.SetMarkerStyle(20);
  hist.SetMarkerSize(0.8);
  hist.GetXaxis()->SetTitle("hit timing [1/6 BCO]");
  hist.GetYaxis()->SetTitle("Number of INTT hits (event averaged)");

  const double y_max = std::max(1.0, hist.GetMaximum() * 1.25);
  hist.GetYaxis()->SetRangeUser(std::min(-5.0, hist.GetMinimum() * 1.2), y_max);

  if (config.draw_full_x_range) {
    hist.GetXaxis()->SetRangeUser(hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
  } else {
    hist.GetXaxis()->SetRangeUser(hist.GetXaxis()->GetBinLowEdge(config.origin_bin),
                                  hist.GetXaxis()->GetBinUpEdge(config.fit_last_bin));
  }

  hist.Draw("E1");
  if (draw_function != nullptr) {
    draw_function->Draw("same");
  }
  if (fine_model_hist != nullptr) {
    fine_model_hist->Draw("hist same");
  }
  DrawFitGuides(hist, config, y_max);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.DrawLatex(0.16, 0.86, multiplicity.label);

  if (summary.status == 0) {
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.16, 0.78,
                    Form("A = %.3f #pm %.3f", summary.amplitude, summary.amplitude_error));
    latex.DrawLatex(0.16, 0.70,
                    Form("#lambda = %.3f #pm %.3f", summary.lambda, summary.lambda_error));
    latex.DrawLatex(0.16, 0.62,
                    Form("#chi^{2}/ndf = %.2f / %d", summary.chi2, summary.ndf));
  } else {
    latex.SetTextColor(kRed + 1);
    latex.DrawLatex(0.16, 0.74, Form("Fit failed (status = %d)", summary.status));
    latex.SetTextColor(kBlack);
  }

  if (draw_legend) {
    TLegend legend(0.58, 0.72, 0.88, 0.90);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&hist, "h_solved_n", "lep");
    if (draw_function != nullptr) {
      legend.AddEntry(draw_function, "Coarse-bin fit", "l");
    }
    if (fine_model_hist != nullptr) {
      legend.AddEntry(fine_model_hist, "Underlying fine model", "l");
    }
    legend.Draw();
  }
}

void WriteSummary(const FitConfig& config, const std::vector<FitSummary>& summaries) {
  const TString output_path = Form("%s/solved_n_poisson_subbin_fit_%s.txt",
                                   config.output_dir.c_str(),
                                   config.group_label.c_str());
  std::ofstream output(output_path.Data());
  if (!output) {
    std::cerr << "ERROR: Could not write " << output_path << std::endl;
    return;
  }

  output << "group_label " << config.group_label << '\n';
  output << "origin_bin " << config.origin_bin << '\n';
  output << "fit_first_bin " << config.fit_first_bin << '\n';
  output << "fit_last_bin " << config.fit_last_bin << '\n';
  output << "subbin_factor " << config.subbin_factor << '\n';
  output << '\n';

  output << std::fixed << std::setprecision(6);
  for (const FitSummary& summary : summaries) {
    output << summary.label << '\n';
    output << "  status " << summary.status << '\n';
    output << "  amplitude " << summary.amplitude << " " << summary.amplitude_error << '\n';
    output << "  lambda " << summary.lambda << " " << summary.lambda_error << '\n';
    output << "  chi2_ndf " << summary.chi2 << " " << summary.ndf << '\n';
  }

  std::cout << "Saved " << output_path << std::endl;
}

}  // namespace

void FitSolvedNPeakPoissonSubbin(const char* group_label = "Group2",
                                 const char* base_output_dir =
                                     "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output",
                                 const char* output_dir =
                                     "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output") {
  SetsPhenixStyle();
  gROOT->SetBatch(kTRUE);

  FitConfig config;
  config.group_label = group_label;
  config.base_output_dir = base_output_dir;
  config.output_dir = output_dir;

  gSystem->mkdir(config.output_dir.c_str(), true);

  const std::vector<MultiplicityConfig> multiplicities = BuildMultiplicityConfigs();
  std::vector<FitSummary> summaries;
  std::vector<std::unique_ptr<TH1D>> histograms;
  std::vector<std::unique_ptr<TF1>> draw_functions;
  std::vector<std::unique_ptr<TH1D>> fine_model_histograms;
  summaries.reserve(multiplicities.size());

  for (const MultiplicityConfig& multiplicity : multiplicities) {
    std::unique_ptr<TH1D> hist(LoadSolvedHistogram(config, multiplicity));
    if (!hist) {
      return;
    }
    if (!ValidateConfig(config, *hist)) {
      return;
    }

    TF1* draw_function = nullptr;
    TH1D* fine_model_hist = nullptr;
    summaries.push_back(
        FitSingleHistogram(*hist, config, multiplicity, draw_function, fine_model_hist));

    histograms.push_back(std::move(hist));
    draw_functions.emplace_back(draw_function);
    fine_model_histograms.emplace_back(fine_model_hist);
  }

  TCanvas canvas(Form("c_poisson_subbin_fit_%s", config.group_label.c_str()),
                 config.group_label.c_str(),
                 900,
                 1100);
  canvas.Divide(1, 3);

  for (std::size_t index = 0; index < histograms.size(); ++index) {
    canvas.cd(static_cast<int>(index) + 1);
    gPad->SetMargin(0.12, 0.04, 0.14, 0.08);
    DrawSinglePad(*histograms[index],
                  config,
                  multiplicities[index],
                  summaries[index],
                  draw_functions[index].get(),
                  fine_model_histograms[index].get(),
                  index == 0);
    if (index == 0) {
      PutCaption(-0.1, 0.65, 0.94);
    }
  }

  const TString output_pdf = Form("%s/solved_n_poisson_subbin_fit_%s.pdf",
                                  config.output_dir.c_str(),
                                  config.group_label.c_str());
  canvas.Print(output_pdf);
  std::cout << "Saved " << output_pdf << std::endl;

  WriteSummary(config, summaries);
}

void PutCaption(double first_margin, double pos_x, double pos_y) {
  TLatex tex;
  const double line_height = 0.055;

  pos_y -= line_height;
  tex.DrawLatexNDC(pos_x, pos_y, "#it{05/13/2026}");
  pos_y -= line_height;
  tex.DrawLatexNDC(pos_x, pos_y, "#it{#bf{sPHENIX}} Internal");
  pos_y -= line_height;
  tex.DrawLatexNDC(pos_x, pos_y, "Run-24 #it{p+p}");
  pos_y -= line_height;
  tex.DrawLatexNDC(pos_x, pos_y, "#sqrt{s} = 200 GeV");
}
