#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dirent.h>

#include "TFile.h"
#include "TH1D.h"
#include "TParameter.h"
#include "TString.h"

namespace {

struct PlateauStats {
  double mean = 0.0;
  double stddev = 0.0;
  int count = 0;
};

bool IsTargetFile(const char* name) {
  const TString filename(name);
  return filename.BeginsWith("run") && filename.EndsWith(".root");
}

PlateauStats CalculatePlateauStats(const TH1D* hist,
                                   int parity,
                                   int peak_bin,
                                   int x_bin_min,
                                   int x_bin_max) {
  std::vector<double> values;

  for (int bin = x_bin_min; bin <= x_bin_max; ++bin) {
    const double x = hist->GetBinCenter(bin);
    if (!(0.0 < x && x < 80.0)) {
      continue;
    }
    if (std::abs(bin - peak_bin) <= 1) {
      continue;
    }
    if ((bin % 2) != parity) {
      continue;
    }
    values.push_back(hist->GetBinContent(bin));
  }

  PlateauStats stats;
  stats.count = values.size();
  if (values.empty()) {
    return stats;
  }

  double sum = 0.0;
  for (double value : values) {
    sum += value;
  }
  stats.mean = sum / values.size();

  double squared_sum = 0.0;
  for (double value : values) {
    const double diff = value - stats.mean;
    squared_sum += diff * diff;
  }
  stats.stddev = std::sqrt(squared_sum / values.size());
  return stats;
}

bool ProcessFile(const std::string& filepath) {
  TFile* file = TFile::Open(filepath.c_str(), "UPDATE");
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << filepath << std::endl;
    delete file;
    return false;
  }

  TH1D* hist = dynamic_cast<TH1D*>(file->Get("h_bco_diff_shifted"));
  if (!hist) {
    std::cerr << "WARNING: h_bco_diff_shifted not found in " << filepath << std::endl;
    file->Close();
    delete file;
    return false;
  }

  const int x_bin_min = hist->FindFixBin(0.0 + 1e-9);
  const int x_bin_max = hist->FindFixBin(80.0 - 1e-9);

  int peak_bin = x_bin_min;
  double peak_value = hist->GetBinContent(x_bin_min);
  for (int bin = x_bin_min; bin <= x_bin_max; ++bin) {
    const double x = hist->GetBinCenter(bin);
    if (!(0.0 < x && x < 80.0)) {
      continue;
    }
    if (hist->GetBinContent(bin) > peak_value) {
      peak_value = hist->GetBinContent(bin);
      peak_bin = bin;
    }
  }

  const PlateauStats even_stats = CalculatePlateauStats(hist, 0, peak_bin, x_bin_min, x_bin_max);
  const PlateauStats odd_stats = CalculatePlateauStats(hist, 1, peak_bin, x_bin_min, x_bin_max);

  if (even_stats.count == 0 || odd_stats.count == 0) {
    std::cerr << "ERROR: Could not estimate both plateau levels in " << filepath << std::endl;
    file->Close();
    delete file;
    return false;
  }

  TH1D* plateau_hist = dynamic_cast<TH1D*>(hist->Clone("h_bco_diff_shifted_plateau"));
  TH1D* subtracted_hist = dynamic_cast<TH1D*>(hist->Clone("h_bco_diff_shifted_minus_plateau"));
  if (!plateau_hist || !subtracted_hist) {
    std::cerr << "ERROR: Failed to clone histograms in " << filepath << std::endl;
    file->Close();
    delete file;
    return false;
  }

  plateau_hist->Reset("ICES");
  plateau_hist->SetTitle("Odd/even plateau estimated from h_bco_diff_shifted");
  subtracted_hist->SetTitle("h_bco_diff_shifted with odd/even plateau subtracted");

  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const bool is_even_bin = (bin % 2) == 0;
    const double plateau_value = is_even_bin ? even_stats.mean : odd_stats.mean;
    const double plateau_error = is_even_bin ? even_stats.stddev : odd_stats.stddev;

    plateau_hist->SetBinContent(bin, plateau_value);
    plateau_hist->SetBinError(bin, plateau_error);

    const double original_value = hist->GetBinContent(bin);
    const double original_error = hist->GetBinError(bin);
    subtracted_hist->SetBinContent(bin, original_value - plateau_value);
    subtracted_hist->SetBinError(
        bin, std::sqrt(original_error * original_error + plateau_error * plateau_error));
  }

  const int peak_low_bin = std::max(1, peak_bin - 1);
  const int peak_high_bin = std::min(hist->GetNbinsX(), peak_bin + 1);

  TParameter<double> even_mean("h_bco_diff_shifted_plateau_even_mean", even_stats.mean);
  TParameter<double> odd_mean("h_bco_diff_shifted_plateau_odd_mean", odd_stats.mean);
  TParameter<double> even_stddev("h_bco_diff_shifted_plateau_even_stddev", even_stats.stddev);
  TParameter<double> odd_stddev("h_bco_diff_shifted_plateau_odd_stddev", odd_stats.stddev);
  TParameter<int> even_count("h_bco_diff_shifted_plateau_even_count", even_stats.count);
  TParameter<int> odd_count("h_bco_diff_shifted_plateau_odd_count", odd_stats.count);
  TParameter<int> peak_bin_param("h_bco_diff_shifted_peak_bin", peak_bin);
  TParameter<int> peak_low_bin_param("h_bco_diff_shifted_peak_excluded_low_bin", peak_low_bin);
  TParameter<int> peak_high_bin_param("h_bco_diff_shifted_peak_excluded_high_bin", peak_high_bin);

  file->cd();
  plateau_hist->Write("", TObject::kOverwrite);
  subtracted_hist->Write("", TObject::kOverwrite);
  even_mean.Write("", TObject::kOverwrite);
  odd_mean.Write("", TObject::kOverwrite);
  even_stddev.Write("", TObject::kOverwrite);
  odd_stddev.Write("", TObject::kOverwrite);
  even_count.Write("", TObject::kOverwrite);
  odd_count.Write("", TObject::kOverwrite);
  peak_bin_param.Write("", TObject::kOverwrite);
  peak_low_bin_param.Write("", TObject::kOverwrite);
  peak_high_bin_param.Write("", TObject::kOverwrite);

  std::cout << filepath
            << " peak_bin=" << peak_bin
            << " even_mean=" << even_stats.mean
            << " even_stddev=" << even_stats.stddev
            << " odd_mean=" << odd_stats.mean
            << " odd_stddev=" << odd_stats.stddev
            << std::endl;

  delete plateau_hist;
  delete subtracted_hist;
  file->Close();
  delete file;
  return true;
}

}  // namespace

void SubtractOddEvenPlateau(
    const char* input_dir =
        "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input") {
  DIR* dir = opendir(input_dir);
  if (!dir) {
    std::cerr << "ERROR: Could not open directory " << input_dir << std::endl;
    return;
  }

  std::vector<std::string> filepaths;
  while (dirent* entry = readdir(dir)) {
    if (!IsTargetFile(entry->d_name)) {
      continue;
    }
    filepaths.push_back(std::string(input_dir) + "/" + entry->d_name);
  }
  closedir(dir);

  std::sort(filepaths.begin(), filepaths.end());

  if (filepaths.empty()) {
    std::cerr << "ERROR: No run*.root files found in " << input_dir << std::endl;
    return;
  }

  int num_processed = 0;
  for (const std::string& filepath : filepaths) {
    if (ProcessFile(filepath)) {
      ++num_processed;
    }
  }

  std::cout << "Processed " << num_processed << " / " << filepaths.size()
            << " files." << std::endl;
}
