#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"

namespace {

const char* kDefaultInputRootPath =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.root";
const char* kDefaultOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";

double ComputeInitialSigma(const TH1D* histogram, int first_bin, int second_bin) {
  double weight_sum = 0.0;
  double mean_sum = 0.0;

  for (int bin = first_bin; bin <= second_bin; ++bin) {
    const double weight = std::max(0.0, histogram->GetBinContent(bin));
    weight_sum += weight;
    mean_sum += weight * histogram->GetBinCenter(bin);
  }

  if (weight_sum <= 0.0) {
    return 1.0;
  }

  const double mean = mean_sum / weight_sum;
  double variance_sum = 0.0;
  for (int bin = first_bin; bin <= second_bin; ++bin) {
    const double weight = std::max(0.0, histogram->GetBinContent(bin));
    const double dx = histogram->GetBinCenter(bin) - mean;
    variance_sum += weight * dx * dx;
  }

  const double sigma = std::sqrt(variance_sum / weight_sum);
  return sigma > 0.0 ? sigma : 1.0;
}

}  // namespace

void FitSolvedNPeak(int first_bin,
                    int second_bin,
                    const char* input_root_path = kDefaultInputRootPath,
                    const char* output_pdf_path = "") {
  if (first_bin > second_bin) {
    std::swap(first_bin, second_bin);
  }

  TFile input_file(input_root_path, "READ");
  if (input_file.IsZombie()) {
    std::cerr << "ERROR: Could not open " << input_root_path << std::endl;
    return;
  }

  TH1D* input_histogram = dynamic_cast<TH1D*>(input_file.Get("h_solved_n"));
  if (input_histogram == nullptr) {
    std::cerr << "ERROR: h_solved_n is missing in " << input_root_path << std::endl;
    return;
  }

  std::unique_ptr<TH1D> histogram(
      static_cast<TH1D*>(input_histogram->Clone("h_solved_n_peak_fit")));
  histogram->SetDirectory(nullptr);
  input_file.Close();

  const int num_bins = histogram->GetNbinsX();
  if (first_bin < 1 || second_bin > num_bins) {
    std::cerr << "ERROR: Requested fit range bins [" << first_bin << ", "
              << second_bin << "] is outside histogram bins [1, " << num_bins
              << "]." << std::endl;
    return;
  }

  int peak_bin = first_bin;
  double peak_value = histogram->GetBinContent(first_bin);
  for (int bin = first_bin + 1; bin <= second_bin; ++bin) {
    if (histogram->GetBinContent(bin) > peak_value) {
      peak_value = histogram->GetBinContent(bin);
      peak_bin = bin;
    }
  }

  const double fit_x_min = histogram->GetXaxis()->GetBinLowEdge(first_bin);
  const double fit_x_max = histogram->GetXaxis()->GetBinUpEdge(second_bin);
  const double initial_amplitude = std::max(peak_value, 1e-6);
  const double initial_mean = histogram->GetBinCenter(peak_bin);
  const double initial_sigma = ComputeInitialSigma(histogram.get(), first_bin, second_bin);

  TF1 fit_function("fit_gaus", "gaus", fit_x_min, fit_x_max);
  fit_function.SetParameters(initial_amplitude, initial_mean, initial_sigma);
  fit_function.SetLineColor(kRed + 1);
  fit_function.SetLineWidth(2);

  TFitResultPtr fit_result = histogram->Fit(&fit_function, "RSQ0");
  const int fit_status = fit_result.Get() ? fit_result->Status() : -1;
  if (!fit_result.Get() || fit_status != 0) {
    std::cerr << "ERROR: Gaussian fit failed for bins [" << first_bin << ", "
              << second_bin << "]. Fit status = " << fit_status << std::endl;
    return;
  }

  const double sigma = std::abs(fit_function.GetParameter(2));
  const double sigma_error = fit_function.GetParError(2);

  gSystem->mkdir(kDefaultOutputDir, true);

  TString output_path = output_pdf_path;
  if (output_path.IsNull()) {
    output_path = Form(
        "%s/solved_n_gaus_fit_bin%dto%d.pdf", kDefaultOutputDir, first_bin, second_bin);
  }
  gSystem->mkdir(gSystem->DirName(output_path), true);

  TCanvas canvas("c_solved_n_peak_fit", "c_solved_n_peak_fit", 900, 700);
  canvas.SetMargin(0.12, 0.04, 0.12, 0.08);

  histogram->SetStats(0);
  histogram->SetLineColor(kBlack);
  histogram->SetLineWidth(2);
  histogram->SetMarkerStyle(20);
  histogram->SetMarkerSize(0.9);
  histogram->SetTitle("h_solved_n Gaussian Peak Fit;bin index i;solved_n[i]");
  histogram->Draw("E1");
  fit_function.Draw("SAME");

  const double y_max = histogram->GetMaximum();
  TLine left_line(first_bin - 0.5, 0.0, first_bin - 0.5, y_max * 1.05);
  TLine right_line(second_bin + 0.5, 0.0, second_bin + 0.5, y_max * 1.05);
  left_line.SetLineStyle(2);
  right_line.SetLineStyle(2);
  left_line.SetLineColor(kBlue + 2);
  right_line.SetLineColor(kBlue + 2);
  left_line.Draw();
  right_line.Draw();

  TLatex caption;
  caption.SetNDC();
  caption.SetTextSize(0.035);
  caption.DrawLatex(
      0.13, 0.92,
      Form("Gaussian fit on bins %d-%d: #sigma = %.4f #pm %.4f", first_bin,
           second_bin, sigma, sigma_error));

  canvas.SaveAs(output_path);

  std::cout << "Fit range: bins [" << first_bin << ", " << second_bin << "]"
            << std::endl;
  std::cout << "Sigma = " << sigma << " +/- " << sigma_error << std::endl;
  std::cout << "Saved PDF: " << output_path << std::endl;
}
