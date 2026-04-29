#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

struct ObservationIndex {
  int relative_bin_position;
  int hist_bin;
  int delay;
  int run;
};

struct OffsetTerm {
  std::string label;
  double coefficient;
};

struct EquationSpec {
  ObservationIndex observation;
  std::vector<int> unknown_indices;
  std::vector<OffsetTerm> offset_terms;
};

struct SolvedParameters {
  std::vector<double> fine_values;
  std::map<std::string, double> offset_values;
};

struct Contribution {
  std::string label;
  double value;
  int color;
};

namespace {

// const char* kDefaultEquationSpecPath = "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/macro/LinearEquationsToBeSolved.c";
const char* kDefaultEquationSpecPath = "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/macro/LinearEquationsToBeSolved_5BCO_18FineBins_scan1.c";
const char* kDefaultSolvedRootPath =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.root";
const char* kDefaultInputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input";
const char* kDefaultOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/shifted_sum_reconstruction";

const int kBeforePeakEvenColor = kOrange + 1;
const int kBeforePeakOddColor = kAzure + 1;
const int kRegularOffsetColor = kGreen + 2;
const int kBoundaryLineColor = kGray + 2;

int GetRunFromDelayValue(int delay) {
  std::vector<std::pair<int, int> > run_L1delay_list_scan6 = {
      {43291, 127}, {43288, 126}, {43285, 125}, {43283, 124}, {43282, 123},
      {43280, 122}, {43278, 121}, 
      {43276, 120}, {43313, 119}
    };
  std::vector<std::pair<int, int> > run_L1delay_list_scan7 = {
      {43408, 120}, {43410, 119}, {43412, 118}, {43413, 117}, {43414, 116},
      {43415, 115}, {43417, 114}, {43421, 113}, {43426, 112}, {43441, 110},
      {43434, 109}, {43436, 108}, {43438, 107}, {43440, 106}};

  for (int i = 0; i < static_cast<int>(run_L1delay_list_scan6.size()); ++i) {
    if (delay == run_L1delay_list_scan6[i].second) {
      return run_L1delay_list_scan6[i].first;
    }
  }

  for (int i = 0; i < static_cast<int>(run_L1delay_list_scan7.size()); ++i) {
    if (delay == run_L1delay_list_scan7[i].second) {
      return run_L1delay_list_scan7[i].first;
    }
  }
  return -1;
}

int GetDelayFromRunValue(int run) {
  std::vector<std::pair<int, int> > run_L1delay_list_scan6 = {
      {43291, 127}, {43288, 126}, {43285, 125}, {43283, 124}, {43282, 123},
      {43280, 122}, {43278, 121}, {43276, 120}, {43313, 119}};
  std::vector<std::pair<int, int> > run_L1delay_list_scan7 = {
      {43408, 120}, {43410, 119}, {43412, 118}, {43413, 117}, {43414, 116},
      {43415, 115}, {43417, 114}, {43421, 113}, {43426, 112}, {43441, 110},
      {43434, 109}, {43436, 108}, {43438, 107}, {43440, 106}};

  for (int i = 0; i < static_cast<int>(run_L1delay_list_scan6.size()); ++i) {
    if (run == run_L1delay_list_scan6[i].first) {
      return run_L1delay_list_scan6[i].second;
    }
  }

  for (int i = 0; i < static_cast<int>(run_L1delay_list_scan7.size()); ++i) {
    if (run == run_L1delay_list_scan7[i].first) {
      return run_L1delay_list_scan7[i].second;
    }
  }
  return -1;
}

std::vector<EquationSpec> LoadEquationSpecs(const std::string& path) {
  std::ifstream input(path.c_str());
  if (!input) {
    std::cerr << "ERROR: Could not open equation spec file " << path << std::endl;
    return {};
  }

  const std::regex prefix_regex(
      R"(^N_\{(?:relative_bin_position|window_bin)=(-?\d+), hist_bin=(\d+), delay=(\d+)\} = (.*)$)");
  const std::regex unknown_regex(R"(n_(\d+))");
  const std::regex offset_term_regex(
      R"(((?:[+-]?\d+(?:\.\d+)?)\s*\*\s*)?(offset_[A-Za-z0-9_]+))");

  std::vector<EquationSpec> equation_specs;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.find("No constraints here") != std::string::npos) {
      continue;
    }

    std::smatch prefix_match;
    if (!std::regex_match(line, prefix_match, prefix_regex)) {
      std::cerr << "WARNING: Could not parse equation line: " << line << std::endl;
      continue;
    }

    EquationSpec spec;
    spec.observation.relative_bin_position = std::stoi(prefix_match[1].str());
    spec.observation.hist_bin = std::stoi(prefix_match[2].str());
    spec.observation.delay = std::stoi(prefix_match[3].str());
    spec.observation.run = GetRunFromDelayValue(spec.observation.delay);

    std::string right_hand_side = prefix_match[4].str();
    const std::size_t value_separator = right_hand_side.rfind(" = ");
    if (value_separator != std::string::npos) {
      right_hand_side = right_hand_side.substr(0, value_separator);
    }

    for (std::sregex_iterator it(right_hand_side.begin(), right_hand_side.end(), unknown_regex), end;
         it != end; ++it) {
      spec.unknown_indices.push_back(std::stoi((*it)[1].str()) - 1);
    }

    for (std::sregex_iterator it(right_hand_side.begin(), right_hand_side.end(), offset_term_regex), end;
         it != end; ++it) {
      std::string coefficient_text = (*it)[1].matched ? (*it)[1].str() : "";
      std::string label = (*it)[2].str();
      double coefficient = 1.0;
      if (!coefficient_text.empty()) {
        coefficient_text.erase(
            std::remove(coefficient_text.begin(), coefficient_text.end(), ' '),
            coefficient_text.end());
        if (!coefficient_text.empty() && coefficient_text.back() == '*') {
          coefficient_text.pop_back();
        }
        coefficient = std::stod(coefficient_text);
      }
      spec.offset_terms.push_back({label, coefficient});
    }

    std::sort(spec.unknown_indices.begin(), spec.unknown_indices.end());
    equation_specs.push_back(spec);
  }

  std::sort(equation_specs.begin(), equation_specs.end(),
            [](const EquationSpec& lhs, const EquationSpec& rhs) {
              if (lhs.observation.run != rhs.observation.run) {
                return lhs.observation.run < rhs.observation.run;
              }
              if (lhs.observation.hist_bin != rhs.observation.hist_bin) {
                return lhs.observation.hist_bin < rhs.observation.hist_bin;
              }
              return lhs.observation.relative_bin_position <
                     rhs.observation.relative_bin_position;
            });
  return equation_specs;
}

SolvedParameters LoadSolvedParameters(const std::string& solved_root_path) {
  SolvedParameters solved_parameters;

  TFile input_file(solved_root_path.c_str(), "READ");
  if (input_file.IsZombie()) {
    std::cerr << "ERROR: Could not open " << solved_root_path << std::endl;
    return solved_parameters;
  }

  TH1D* h_solved_n = dynamic_cast<TH1D*>(input_file.Get("h_solved_n"));
  if (h_solved_n == nullptr) {
    std::cerr << "ERROR: h_solved_n is missing in " << solved_root_path << std::endl;
    return solved_parameters;
  }

  solved_parameters.fine_values.resize(h_solved_n->GetNbinsX(), 0.0);
  for (int bin = 1; bin <= h_solved_n->GetNbinsX(); ++bin) {
    solved_parameters.fine_values[bin - 1] = h_solved_n->GetBinContent(bin);
  }

  TH1D* h_solved_offsets = dynamic_cast<TH1D*>(input_file.Get("h_solved_offsets"));
  if (h_solved_offsets == nullptr) {
    std::cerr << "ERROR: h_solved_offsets is missing in " << solved_root_path << std::endl;
    std::cerr << "Run SolveShiftedSum.C again after updating it so offsets are saved."
              << std::endl;
    solved_parameters.fine_values.clear();
    return solved_parameters;
  }

  for (int bin = 1; bin <= h_solved_offsets->GetNbinsX(); ++bin) {
    solved_parameters.offset_values[h_solved_offsets->GetXaxis()->GetBinLabel(bin)] =
        h_solved_offsets->GetBinContent(bin);
  }

  return solved_parameters;
}

std::map<int, std::vector<EquationSpec> > GroupEquationsByRun(
    const std::vector<EquationSpec>& equation_specs) {
  std::map<int, std::vector<EquationSpec> > grouped;
  for (int i = 0; i < static_cast<int>(equation_specs.size()); ++i) {
    if (equation_specs[i].observation.run < 0) {
      std::cerr << "WARNING: Unknown run for delay "
                << equation_specs[i].observation.delay << std::endl;
      continue;
    }
    grouped[equation_specs[i].observation.run].push_back(equation_specs[i]);
  }
  return grouped;
}

int GetFineBinColor(int fine_bin_index) {
  static const int colors[] = {
      kRed + 1,    kMagenta + 1, kViolet - 1, kBlue + 1,
      kCyan + 1,   kTeal + 1,    kSpring + 5, kPink + 9};
  return colors[fine_bin_index % (sizeof(colors) / sizeof(colors[0]))];
}

std::string GetFineBinLabel(int fine_bin_index) {
  return Form("h_solved_n bin %d", fine_bin_index + 1);
}

int GetOffsetColor(const std::string& label) {
  if (label.find("before_peak_even") != std::string::npos) {
    return kBeforePeakEvenColor;
  }
  if (label.find("before_peak_odd") != std::string::npos) {
    return kBeforePeakOddColor;
  }
  return kRegularOffsetColor;
}

std::string GetOffsetLegendLabel(const std::string& label) {
  if (label.find("before_peak_even") != std::string::npos) {
    return "before_peak_even";
  }
  if (label.find("before_peak_odd") != std::string::npos) {
    return "before_peak_odd";
  }
  return "offset";
}

TH1D* LoadRunHistogram(const std::string& input_dir, int run) {
  const std::string path = input_dir + Form("/run%d.root", run);
  TFile* input_file = TFile::Open(path.c_str(), "READ");
  if (input_file == nullptr || input_file->IsZombie()) {
    std::cerr << "WARNING: Could not open " << path << std::endl;
    delete input_file;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(input_file->Get("h_bco_diff_shifted"));
  if (hist == nullptr) {
    std::cerr << "WARNING: h_bco_diff_shifted is missing in " << path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  TH1D* cloned_hist = dynamic_cast<TH1D*>(hist->Clone(Form("h_bco_diff_shifted_run%d", run)));
  cloned_hist->SetDirectory(nullptr);
  input_file->Close();
  delete input_file;
  return cloned_hist;
}

double ComputeMaximumMagnitude(const TH1D* hist) {
  double max_magnitude = 0.0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    max_magnitude = std::max(max_magnitude, std::abs(hist->GetBinContent(bin)));
  }
  return max_magnitude;
}

double ComputeMinimumPositive(const TH1D* hist) {
  double min_positive = -1.0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double value = hist->GetBinContent(bin);
    if (value <= 0.0) {
      continue;
    }
    if (min_positive < 0.0 || value < min_positive) {
      min_positive = value;
    }
  }
  return min_positive;
}

void DrawRunDelayLabel(int run) {
  const int delay = GetDelayFromRunValue(run);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(13);
  latex.SetTextSize(0.032);
  if (delay >= 0) {
    latex.DrawLatex(0.11, 0.965, Form("run %d, L1delay %d", run, delay));
  } else {
    latex.DrawLatex(0.11, 0.965, Form("run %d, L1delay unknown", run));
  }
}

void DrawContributionArrowSequence(double x_low,
                                   double x_high,
                                   const std::vector<Contribution>& contributions,
                                   double min_draw_y) {
  if (contributions.empty()) {
    return;
  }

  const double width = x_high - x_low;
  const double margin = 0.14 * width;
  const double usable_x_low = x_low + margin;
  const double usable_x_high = x_high - margin;
  const int num_steps = static_cast<int>(contributions.size());
  const double delta_x =
      (num_steps > 1) ? (usable_x_high - usable_x_low) / (num_steps - 1) : 0.0;

  double cumulative = 0.0;
  double previous_x = usable_x_low;

  for (int i = 0; i < num_steps; ++i) {
    const Contribution& contribution = contributions[i];
    const double current_x = usable_x_low + delta_x * i;
    const double next_cumulative = cumulative + contribution.value;
    const double draw_start = std::max(cumulative, min_draw_y);
    const double draw_end = std::max(next_cumulative, min_draw_y);

    if (i > 0) {
      TLine* connector = new TLine(previous_x, draw_start, current_x, draw_start);
      connector->SetLineColor(kGray + 1);
      connector->SetLineStyle(1);
      connector->SetLineWidth(1);
      connector->Draw();
    }

    if (draw_end > draw_start) {
      TArrow* arrow = new TArrow(current_x, draw_start, current_x, draw_end, 0.015, "|>");
      arrow->SetLineColor(contribution.color);
      arrow->SetFillColor(contribution.color);
      arrow->SetLineWidth(2);
      arrow->Draw();
    } else if (draw_end < draw_start) {
      TArrow* arrow = new TArrow(current_x, draw_start, current_x, draw_end, 0.015, "|>");
      arrow->SetLineColor(contribution.color);
      arrow->SetFillColor(contribution.color);
      arrow->SetLineWidth(2);
      arrow->Draw();
    } else {
      TLine* clipped = new TLine(current_x - 0.015 * width, draw_start,
                                 current_x + 0.015 * width, draw_start);
      clipped->SetLineColor(contribution.color);
      clipped->SetLineWidth(2);
      clipped->Draw();
    }

    cumulative = next_cumulative;
    previous_x = current_x;
  }

  const double draw_final = std::max(cumulative, min_draw_y);
  TLine* final_level = new TLine(usable_x_high, draw_final, x_high - 0.04 * width, draw_final);
  final_level->SetLineColor(kGray + 1);
  final_level->SetLineWidth(1);
  final_level->Draw();
}

void DrawLegend(std::vector<TH1D*>* dummies,
                const std::set<int>& used_fine_bins,
                bool has_before_peak_even,
                bool has_before_peak_odd,
                bool has_regular_offset,
                bool include_boundary_entry) {
  TLegend* legend = new TLegend(0.67, 0.60, 0.93, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.030);

  if (has_before_peak_even) {
    TH1D* dummy = new TH1D(Form("dummy_before_even_%zu", dummies->size()), "", 1, 0.0, 1.0);
    dummy->SetLineWidth(3);
    dummy->SetLineColor(kBeforePeakEvenColor);
    dummies->push_back(dummy);
    legend->AddEntry(dummy, "before_peak_even", "l");
  }

  if (has_before_peak_odd) {
    TH1D* dummy = new TH1D(Form("dummy_before_odd_%zu", dummies->size()), "", 1, 0.0, 1.0);
    dummy->SetLineWidth(3);
    dummy->SetLineColor(kBeforePeakOddColor);
    dummies->push_back(dummy);
    legend->AddEntry(dummy, "before_peak_odd", "l");
  }

  if (has_regular_offset) {
    TH1D* dummy = new TH1D(Form("dummy_offset_%zu", dummies->size()), "", 1, 0.0, 1.0);
    dummy->SetLineWidth(3);
    dummy->SetLineColor(kRegularOffsetColor);
    dummies->push_back(dummy);
    legend->AddEntry(dummy, "offset", "l");
  }

  for (std::set<int>::const_iterator it = used_fine_bins.begin(); it != used_fine_bins.end();
       ++it) {
    const int fine_bin_index = *it;
    TH1D* dummy = new TH1D(Form("dummy_fine_bin_%d_%zu", fine_bin_index, dummies->size()), "",
                           1, 0.0, 1.0);
    dummy->SetLineWidth(3);
    dummy->SetLineColor(GetFineBinColor(fine_bin_index));
    dummies->push_back(dummy);
    legend->AddEntry(dummy, GetFineBinLabel(fine_bin_index).c_str(), "l");
  }

  if (include_boundary_entry) {
    TH1D* dummy = new TH1D(Form("dummy_boundary_%zu", dummies->size()), "", 1, 0.0, 1.0);
    dummy->SetLineColor(kBoundaryLineColor);
    dummy->SetLineStyle(2);
    dummy->SetLineWidth(2);
    dummies->push_back(dummy);
    legend->AddEntry(dummy, "coarse-bin boundary", "l");
  }

  legend->Draw();
}

void SaveFigure1(const TH1D* h_solved_n,
                 const std::vector<EquationSpec>& run_equations,
                 int run,
                 const std::string& output_dir) {
  TCanvas canvas(Form("c_solved_n_partition_run%d", run), "", 1400, 700);

  const double max_magnitude = std::max(1.0, ComputeMaximumMagnitude(h_solved_n));
  const double y_min = -0.15 * max_magnitude;
  const double y_max = 1.25 * max_magnitude;

  TH1D* frame = dynamic_cast<TH1D*>(h_solved_n->Clone(Form("h_solved_n_frame_run%d", run)));
  frame->Reset("ICES");
  frame->SetStats(0);
  frame->SetTitle(Form("run %d: solved fine bins grouped into coarse bins;fine bin index;solved_n", run));
  frame->SetMinimum(y_min);
  frame->SetMaximum(y_max);
  frame->Draw("AXIS");

  std::set<double> boundaries;
  std::set<int> used_fine_bins;
  for (int i = 0; i < static_cast<int>(run_equations.size()); ++i) {
    const EquationSpec& equation = run_equations[i];
    if (equation.unknown_indices.empty()) {
      continue;
    }

    const int first_hist_bin = equation.unknown_indices.front() + 1;
    const int last_hist_bin = equation.unknown_indices.back() + 1;
    boundaries.insert(h_solved_n->GetXaxis()->GetBinLowEdge(first_hist_bin));
    boundaries.insert(h_solved_n->GetXaxis()->GetBinUpEdge(last_hist_bin));

    for (int rank = 0; rank < static_cast<int>(equation.unknown_indices.size()); ++rank) {
      const int fine_bin_index = equation.unknown_indices[rank];
      const int fine_hist_bin = fine_bin_index + 1;
      const double value = h_solved_n->GetBinContent(fine_hist_bin);
      const double x_low = h_solved_n->GetXaxis()->GetBinLowEdge(fine_hist_bin);
      const double x_high = h_solved_n->GetXaxis()->GetBinUpEdge(fine_hist_bin);
      TBox* box = new TBox(x_low, std::min(0.0, value), x_high, std::max(0.0, value));
      box->SetFillColorAlpha(GetFineBinColor(fine_bin_index), 0.55);
      box->SetLineColor(GetFineBinColor(fine_bin_index));
      box->Draw();
      used_fine_bins.insert(fine_bin_index);
    }
  }

  TH1D* overlay = dynamic_cast<TH1D*>(h_solved_n->Clone(Form("h_solved_n_overlay_run%d", run)));
  overlay->SetStats(0);
  overlay->SetFillStyle(0);
  overlay->SetLineColor(kBlack);
  overlay->SetLineWidth(2);
  overlay->Draw("HIST SAME");

  for (std::set<double>::const_iterator it = boundaries.begin(); it != boundaries.end(); ++it) {
    TLine* line = new TLine(*it, y_min, *it, y_max);
    line->SetLineColor(kBoundaryLineColor);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
  }

  TLatex latex;
  latex.SetTextAlign(22);
  latex.SetTextSize(0.030);
  for (int i = 0; i < static_cast<int>(run_equations.size()); ++i) {
    const EquationSpec& equation = run_equations[i];
    if (equation.unknown_indices.empty()) {
      continue;
    }
    const double x_center =
        0.5 * (equation.unknown_indices.front() + 1 + equation.unknown_indices.back() + 1);
    latex.DrawLatex(x_center, y_max * 0.96, Form("hist bin %d", equation.observation.hist_bin));
  }

  DrawRunDelayLabel(run);

  canvas.SaveAs((output_dir + Form("/run%d_solved_n_partition.pdf", run)).c_str());

  delete frame;
  delete overlay;
}

void SaveFigure2(const TH1D* h_bco_diff_shifted,
                 const std::vector<EquationSpec>& run_equations,
                 const SolvedParameters& solved_parameters,
                 int run,
                 const std::string& output_dir,
                 bool use_log_y) {
  std::map<int, std::vector<Contribution> > contributions_by_hist_bin;
  bool has_before_peak_even = false;
  bool has_before_peak_odd = false;
  bool has_regular_offset = false;
  std::set<int> used_fine_bins;
  double reconstructed_max_magnitude = 0.0;
  double reconstructed_max_positive = 0.0;
  double reconstructed_min_positive = -1.0;

  for (int i = 0; i < static_cast<int>(run_equations.size()); ++i) {
    const EquationSpec& equation = run_equations[i];
    std::vector<Contribution>& contributions =
        contributions_by_hist_bin[equation.observation.hist_bin];

    for (int term_index = 0; term_index < static_cast<int>(equation.offset_terms.size());
         ++term_index) {
      const OffsetTerm& offset_term = equation.offset_terms[term_index];
      const std::map<std::string, double>::const_iterator offset_it =
          solved_parameters.offset_values.find(offset_term.label);
      if (offset_it == solved_parameters.offset_values.end()) {
        std::cerr << "WARNING: Missing solved value for " << offset_term.label
                  << " while processing run " << run << std::endl;
        continue;
      }

      Contribution contribution;
      contribution.label = GetOffsetLegendLabel(offset_term.label);
      contribution.value = offset_term.coefficient * offset_it->second;
      contribution.color = GetOffsetColor(offset_term.label);
      contributions.push_back(contribution);

      if (contribution.label == "before_peak_even") {
        has_before_peak_even = true;
      } else if (contribution.label == "before_peak_odd") {
        has_before_peak_odd = true;
      } else {
        has_regular_offset = true;
      }
    }

    for (int rank = 0; rank < static_cast<int>(equation.unknown_indices.size()); ++rank) {
      const int unknown_index = equation.unknown_indices[rank];
      if (unknown_index < 0 ||
          unknown_index >= static_cast<int>(solved_parameters.fine_values.size())) {
        std::cerr << "WARNING: Fine-bin index out of range while processing run " << run
                  << std::endl;
        continue;
      }

      Contribution contribution;
      contribution.label = GetFineBinLabel(unknown_index);
      contribution.value = solved_parameters.fine_values[unknown_index];
      contribution.color = GetFineBinColor(unknown_index);
      contributions.push_back(contribution);
      used_fine_bins.insert(unknown_index);
    }

    double reconstructed_sum = 0.0;
    double cumulative = 0.0;
    for (int j = 0; j < static_cast<int>(contributions.size()); ++j) {
      reconstructed_sum += contributions[j].value;
      cumulative += contributions[j].value;
      if (cumulative > 0.0) {
        reconstructed_max_positive = std::max(reconstructed_max_positive, cumulative);
        if (reconstructed_min_positive < 0.0 || cumulative < reconstructed_min_positive) {
          reconstructed_min_positive = cumulative;
        }
      }
    }
    reconstructed_max_magnitude =
        std::max(reconstructed_max_magnitude, std::abs(reconstructed_sum));
  }

  TCanvas canvas(Form("c_h_bco_diff_shifted_run%d", run), "", 1400, 700);
  canvas.SetLogy(use_log_y);

  double y_min = 0.0;
  double y_max = 0.0;
  if (use_log_y) {
    y_min = 0.01;
    y_max = 2.0 * 100.0;
  } else {
    y_min = -7.5;
    y_max = 90.0;
  }

  TH1D* frame =
      dynamic_cast<TH1D*>(h_bco_diff_shifted->Clone(Form("h_bco_diff_shifted_frame_run%d", run)));
  frame->Reset("ICES");
  frame->SetStats(0);
  frame->SetTitle(Form("run %d: coarse-bin reconstruction from solved fine bins;h_{bco_diff_shifted} bin;events", run));
  frame->SetMinimum(y_min);
  frame->SetMaximum(y_max);
  frame->GetXaxis()->SetRangeUser(0,7);
  frame->Draw("AXIS");

  for (std::map<int, std::vector<Contribution> >::const_iterator it =
           contributions_by_hist_bin.begin();
       it != contributions_by_hist_bin.end(); ++it) {
    const int hist_bin = it->first;
    if (hist_bin < 1 || hist_bin > h_bco_diff_shifted->GetNbinsX()) {
      continue;
    }

    const double x_low = h_bco_diff_shifted->GetXaxis()->GetBinLowEdge(hist_bin);
    const double x_high = h_bco_diff_shifted->GetXaxis()->GetBinUpEdge(hist_bin);
    DrawContributionArrowSequence(x_low, x_high, it->second,
                                  use_log_y ? y_min : std::numeric_limits<double>::lowest());
  }

  TH1D* overlay = dynamic_cast<TH1D*>(
      h_bco_diff_shifted->Clone(Form("h_bco_diff_shifted_overlay_run%d", run)));
  overlay->SetStats(0);
  overlay->SetFillStyle(0);
  overlay->SetLineColor(kBlack);
  overlay->SetLineWidth(2);
  overlay->SetMarkerStyle(20);
  overlay->SetMarkerSize(0.9);
  overlay->Draw("HIST SAME");
  overlay->Draw("E1 SAME");

  DrawRunDelayLabel(run);

  canvas.SaveAs(
      (output_dir +
       Form("/run%d_h_bco_diff_shifted_reconstruction_%s.pdf", run,
            use_log_y ? "log" : "linear"))
          .c_str());

  delete frame;
  delete overlay;
}

}  // namespace

void VisualizeShiftedSumReconstruction(
    const char* equation_spec_path = kDefaultEquationSpecPath,
    const char* solved_root_path = kDefaultSolvedRootPath,
    const char* input_dir = kDefaultInputDir,
    const char* output_dir = kDefaultOutputDir) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(output_dir, true);

  const std::vector<EquationSpec> equation_specs = LoadEquationSpecs(equation_spec_path);
  if (equation_specs.empty()) {
    std::cerr << "ERROR: No equation specs were loaded." << std::endl;
    return;
  }

  const SolvedParameters solved_parameters = LoadSolvedParameters(solved_root_path);
  if (solved_parameters.fine_values.empty()) {
    std::cerr << "ERROR: Solved parameters are unavailable." << std::endl;
    return;
  }

  TFile solved_file(solved_root_path, "READ");
  TH1D* h_solved_n_in_file = dynamic_cast<TH1D*>(solved_file.Get("h_solved_n"));
  if (h_solved_n_in_file == nullptr) {
    std::cerr << "ERROR: h_solved_n is missing in " << solved_root_path << std::endl;
    return;
  }
  TH1D* h_solved_n = dynamic_cast<TH1D*>(h_solved_n_in_file->Clone("h_solved_n_visualization"));
  h_solved_n->SetDirectory(nullptr);
  solved_file.Close();

  const std::map<int, std::vector<EquationSpec> > equations_by_run =
      GroupEquationsByRun(equation_specs);

  for (std::map<int, std::vector<EquationSpec> >::const_iterator it =
           equations_by_run.begin();
       it != equations_by_run.end(); ++it) {
    const int run = it->first;
    TH1D* h_bco_diff_shifted = LoadRunHistogram(input_dir, run);
    if (h_bco_diff_shifted == nullptr) {
      continue;
    }

    SaveFigure1(h_solved_n, it->second, run, output_dir);
    SaveFigure2(h_bco_diff_shifted, it->second, solved_parameters, run, output_dir, false);
    SaveFigure2(h_bco_diff_shifted, it->second, solved_parameters, run, output_dir, true);
    std::cout << "Created visualization for run " << run << std::endl;

    delete h_bco_diff_shifted;
  }

  delete h_solved_n;
}
