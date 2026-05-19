#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TArrow.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "Utility/sPhenixStyle.C"

void PutCaption(double first_margin = -0.1, double pos_x = 0.2, double pos_y = 0.9);

namespace {

struct MultiplicityConfig {
  const char* dir_suffix;
  const char* legend_label;
  int line_color;
};

const char* kDefaultBaseOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
constexpr double kBcoWidthNs = 106.6;
constexpr double kFineBinWidthNs = kBcoWidthNs / 6.0;
constexpr int kTimeZeroBinLeftEdge = 2;
constexpr int kHighlightedFirstBin = 8;
constexpr int kHighlightedLastBin = 13;
constexpr double kHighlightLineBottom = -5.0;
constexpr double kHighlightLineTop = 160.0;

std::vector<MultiplicityConfig> BuildMultiplicityConfigs() {
  return {
      // {"mul10to100", "10 < trigger timing hits < 100", kRed + 1},
      // {"mul100to200", "100 < trigger timing hits < 200", kBlue + 1},
      // {"mul200to1000", "200 < trigger timing hits < 1000", kGreen + 2},
      {"mul10to100", "Low multiplicity", kRed + 1},
      {"mul100to200", "Middle multiplicity", kBlue + 1},
      {"mul200to1000", "High multiplicity", kGreen + 2},
  };
}

TH1D* LoadSolvedHistogram(const std::string& base_output_dir,
                          const std::string& group_label,
                          const MultiplicityConfig& multiplicity) {
  const TString input_path = Form(
      "%s/shifted_sum_reconstruction_%s_%s/solved_n.root",
      base_output_dir.c_str(),
      group_label.c_str(),
      multiplicity.dir_suffix);

  TFile* input_file = TFile::Open(input_path, "READ");
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << input_path << std::endl;
    delete input_file;
    return nullptr;
  }

  TH1D* input_hist = dynamic_cast<TH1D*>(input_file->Get("h_solved_n"));
  if (input_hist == nullptr) {
    std::cerr << "ERROR: h_solved_n is missing in " << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(
      input_hist->Clone(Form("h_solved_n_%s_%s", group_label.c_str(), multiplicity.dir_suffix)));
  if (hist == nullptr) {
    std::cerr << "ERROR: Failed to clone h_solved_n from " << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  hist->SetDirectory(nullptr);
  input_file->Close();
  delete input_file;
  return hist;
}

TH1D* BuildNanosecondHistogram(const TH1D& source_hist, const TString& histogram_name) {
  const int num_bins = source_hist.GetNbinsX();
  const double x_min =
      (source_hist.GetXaxis()->GetBinLowEdge(1) - kTimeZeroBinLeftEdge + 0.5) * kFineBinWidthNs;
  const double x_max =
      (source_hist.GetXaxis()->GetBinUpEdge(num_bins) - kTimeZeroBinLeftEdge + 0.5) * kFineBinWidthNs;

  TH1D* hist = new TH1D(histogram_name, source_hist.GetTitle(), num_bins, x_min, x_max);
  hist->SetDirectory(nullptr);
  for (int bin = 1; bin <= num_bins; ++bin) {
    hist->SetBinContent(bin, source_hist.GetBinContent(bin));
    hist->SetBinError(bin, source_hist.GetBinError(bin));
  }
  hist->SetEntries(source_hist.GetEntries());
  return hist;
}

void DrawHighlightLines() {
  const double x_min =
      (kHighlightedFirstBin - kTimeZeroBinLeftEdge) * kFineBinWidthNs;
  const double x_max =
      (kHighlightedLastBin - kTimeZeroBinLeftEdge + 1.0) * kFineBinWidthNs;

  TLine* left_line = new TLine(x_min, kHighlightLineBottom, x_min, kHighlightLineTop);
  left_line->SetLineStyle(3);
  left_line->SetLineWidth(2);
  left_line->Draw();

  TLine* right_line = new TLine(x_max, kHighlightLineBottom, x_max, kHighlightLineTop);
  right_line->SetLineStyle(3);
  right_line->SetLineWidth(2);
  right_line->Draw();
}

void DrawHighlightAnnotation() {
  const double x_min =
      (kHighlightedFirstBin - kTimeZeroBinLeftEdge) * kFineBinWidthNs;
  const double x_max =
      (kHighlightedLastBin - kTimeZeroBinLeftEdge + 1.0) * kFineBinWidthNs;
  const double x_center = 0.5 * (x_min + x_max);
  const double x_span = x_max - x_min;
  const double arrow_y = 158.0;
  const double text_y = 175;

  TArrow* left_arrow =
      new TArrow(x_center - 0.10 * x_span, arrow_y, x_min + 0.04 * x_span, arrow_y, 0.018, "|>");
  left_arrow->SetLineWidth(2);
  left_arrow->Draw();

  TArrow* right_arrow =
      new TArrow(x_center + 0.10 * x_span, arrow_y, x_max - 0.04 * x_span, arrow_y, 0.018, "|>");
  right_arrow->SetLineWidth(2);
  right_arrow->Draw();

  TLatex latex;
  latex.SetTextAlign(22);
  // latex.SetTextSize(0.034);
  latex.DrawLatex(x_center, text_y, "1 BCO range at the Run-24 working point");
}

bool SaveStandardPlot(const std::string& group_label,
                      const std::string& base_output_dir,
                      const std::vector<TH1D*>& histograms,
                      const std::vector<MultiplicityConfig>& multiplicities) {
  TCanvas* canvas =
      new TCanvas(Form("c_compare_solved_n_%s", group_label.c_str()), group_label.c_str(), 900, 700);
  TLegend* legend = new TLegend(0.6, 0.7, 0.8, 0.94);

  bool first_hist = true;
  for (std::size_t index = 0; index < histograms.size(); ++index) {
    TH1D* hist = histograms[index];
    hist->SetLineColor(multiplicities[index].line_color);
    hist->GetXaxis()->SetTitle("hit timing [1/6 BCO]");
    hist->GetYaxis()->SetTitle("Number of INTT hits (event averaged)");
    hist->GetYaxis()->SetRangeUser(-5.0, 300.0);

    if (first_hist) {
      hist->Draw("hist");
      first_hist = false;
    } else {
      hist->Draw("hist same");
    }

    legend->AddEntry(hist, multiplicities[index].legend_label, "l");
  }

  PutCaption(-0.1, 0.2, 0.94);
  legend->Draw();

  const TString output_path = Form(
      "%s/fine_hit_timing_distribution_multiplicity_comparison_%s.pdf",
      base_output_dir.c_str(),
      group_label.c_str());
  canvas->Print(output_path);
  std::cout << "Saved " << output_path << std::endl;

  delete legend;
  delete canvas;
  return true;
}

bool SaveNanosecondPlot(const std::string& group_label,
                        const std::string& base_output_dir,
                        const std::vector<TH1D*>& histograms,
                        const std::vector<MultiplicityConfig>& multiplicities) {
  TCanvas* canvas = new TCanvas(Form("c_compare_solved_n_ns_%s", group_label.c_str()),
                                group_label.c_str(),
                                900,
                                700);
  TLegend* legend = new TLegend(0.6, 0.7, 0.8, 0.94);

  std::vector<TH1D*> time_histograms;
  time_histograms.reserve(histograms.size());

  bool first_hist = true;
  for (std::size_t index = 0; index < histograms.size(); ++index) {
    TH1D* time_hist = BuildNanosecondHistogram(
        *histograms[index],
        Form("h_solved_n_ns_%s_%s", group_label.c_str(), multiplicities[index].dir_suffix));
    time_hist->SetLineColor(multiplicities[index].line_color);
    time_hist->GetXaxis()->SetTitle("hit timing [ns]");
    time_hist->GetYaxis()->SetTitle("Number of INTT hits (event averaged)");
    time_hist->GetYaxis()->SetRangeUser(-5.0, 300.0);

    if (first_hist) {
      time_hist->Draw("hist");
      first_hist = false;
    } else {
      time_hist->Draw("hist same");
    }

    legend->AddEntry(time_hist, multiplicities[index].legend_label, "l");
    time_histograms.push_back(time_hist);
  }

  DrawHighlightLines();
  DrawHighlightAnnotation();
  PutCaption(-0.1, 0.2, 0.94);
  legend->Draw();

  const TString output_path = Form(
      "%s/fine_hit_timing_distribution_multiplicity_comparison_ns_%s.pdf",
      base_output_dir.c_str(),
      group_label.c_str());
  canvas->Print(output_path);
  std::cout << "Saved " << output_path << std::endl;

  for (TH1D* hist : time_histograms) {
    delete hist;
  }
  delete legend;
  delete canvas;
  return true;
}

bool DrawGroup(const std::string& group_label,
               const std::string& base_output_dir,
               const std::vector<MultiplicityConfig>& multiplicities) {
  std::vector<TH1D*> histograms;
  histograms.reserve(multiplicities.size());

  for (const MultiplicityConfig& multiplicity : multiplicities) {
    TH1D* hist = LoadSolvedHistogram(base_output_dir, group_label, multiplicity);
    if (hist == nullptr) {
      for (TH1D* loaded_hist : histograms) {
        delete loaded_hist;
      }
      return false;
    }
    histograms.push_back(hist);
  }

  const bool saved_standard =
      SaveStandardPlot(group_label, base_output_dir, histograms, multiplicities);
  const bool saved_ns = SaveNanosecondPlot(group_label, base_output_dir, histograms, multiplicities);

  for (TH1D* hist : histograms) {
    delete hist;
  }
  return saved_standard && saved_ns;
}

}  // namespace

void CompareSolvedNByMultiplicity(const char* base_output_dir = kDefaultBaseOutputDir) {
  SetsPhenixStyle();
  
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(base_output_dir, true);

  const std::vector<MultiplicityConfig> multiplicities = BuildMultiplicityConfigs();
  const std::vector<std::string> groups = {"Group1", "Group1_wo43280", "Group2"};

  int num_success = 0;
  for (const std::string& group_label : groups) {
    if (DrawGroup(group_label, base_output_dir, multiplicities)) {
      ++num_success;
    } else {
      std::cerr << "ERROR: Failed to draw " << group_label << std::endl;
    }
  }

  std::cout << "Created " << num_success << " / " << groups.size()
            << " multiplicity comparison canvases." << std::endl;
}


void PutCaption(double first_margin = -0.1, double pos_x = 0.2, double pos_y = 0.9){
/*from Genki*/
  TLatex* tex = new TLatex();
  double line_height = 0.055; // no need to change
  // double first_margin = -0.1; // the margin from the top of the canvas to the top line (date). You can modify
  // double pos_x = 0.2; // you can modify
  // double pos_y = 1.0  + first_margin; // - line_height;

  /* Date, you may want to change 0.65 to pos_x */
  //tex->DrawLatexNDC( pos_x, pos_y,
  //           string("#it{" + GetDate() + "}").c_str() );

  // sPHENIX Internal or sPHENIX Prelimnary
  pos_y -= line_height;
//   tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Preliminary" );
  // tex->DrawLatexNDC( pos_x, 0.955, "#it{#bf{sPHENIX}} Internal" );
  // tex->DrawLatexNDC( pos_x, 0.955, "#it{#bf{sPHENIX}} Preliminary" );
  //tex->DrawLatexNDC( pos_x, pos_y, string( string("#it{#bf{sPHENIX}} Preliminary ") + "#it{" + GetDate() + "}").c_str() );
  //tex->DrawLatexNDC( pos_x, pos_y, string("#it{#bf{sPHENIX}} Preliminary #it{05/16/2025}").c_str() );

  // tex->DrawLatexNDC( 0.77, 0.955, string("#it{05/13/2026}").c_str() );
  // tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Internal #it{05/13/2026}" );
  tex->DrawLatexNDC( pos_x, pos_y, "#it{05/13/2026}" );
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Internal" );
  
  // p+p 200 GeV
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, ("Run-24 #it{p+p}") );
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, ("#sqrt{s} = 200 GeV") );
  // tex->DrawLatexNDC( pos_x, 0.88, ("#it{p+p} #sqrt{s} = 200 GeV") );


}
