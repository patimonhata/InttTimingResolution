#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TSystem.h"
#include "TVectorD.h"

struct ObservationIndex {
  int relative_bin_position;
  int hist_bin;
  int delay;
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

namespace {
constexpr int kDelayMin = 106;
constexpr int kDelayMax = 127;
const char* kEquationSpecPath =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/macro/LinearEquationsToBeSolved.c";
}

int GetRunFromDelayValue(int delay);
int GetNumUnknowns(const std::vector<EquationSpec>& equation_specs);
std::vector<EquationSpec> LoadEquationSpecs(const std::string& path);
std::vector<std::string> CollectOffsetLabels(const std::vector<EquationSpec>& equation_specs);
void BuildCoefficientMatrix(TMatrixD* matrix_a,
                            int num_unknowns,
                            const std::vector<EquationSpec>& equation_specs,
                            const std::vector<std::string>& offset_labels);
void SaveSolvedHistogram(int num_unknowns,
                         const TVectorD& solved_n,
                         const TMatrixD& covariance_matrix,
                         const std::vector<std::string>& offset_labels);
bool SolveShiftedSumWeighted(const TVectorD& measured_n,
                             const TVectorD& measured_sigma,
                             int num_unknowns,
                             const std::vector<EquationSpec>& equation_specs,
                             const std::vector<std::string>& offset_labels,
                             TVectorD* solved_n,
                             TMatrixD* covariance_matrix,
                             TMatrixD* coefficient_matrix = nullptr);
void PrintSolvedEquations(int num_unknowns,
                          const TMatrixD& matrix_a,
                          const TVectorD& measured_n,
                          const std::vector<EquationSpec>& equation_specs,
                          const std::vector<std::string>& offset_labels);
void PrintSolution(int num_unknowns,
                   const TVectorD& solved_n,
                   const TMatrixD& covariance_matrix,
                   const std::vector<std::string>& offset_labels);

void SolveShiftedSum() {
  const std::vector<EquationSpec> equation_specs = LoadEquationSpecs(kEquationSpecPath);
  const int num_unknowns = GetNumUnknowns(equation_specs);
  const std::vector<std::string> offset_labels = CollectOffsetLabels(equation_specs);
  const int kNumObservations = static_cast<int>(equation_specs.size());

  TVectorD measured_n(kNumObservations);
  TVectorD measured_sigma(kNumObservations);

  measured_n.Zero();
  measured_sigma.Zero();

  for (int row = 0; row < kNumObservations; ++row) {
    measured_sigma(row) = 0.01;
  }

  for (int row = 0; row < kNumObservations; ++row) {
    const ObservationIndex& observation = equation_specs[row].observation;
    const int run = GetRunFromDelayValue(observation.delay);
    std::string filename = Form(
        "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input/run%d.root",
        run);
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "ERROR: Could not open " << filename << std::endl;
      exit(1);
    }

    TH1D* h_bco_diff_shifted = (TH1D*)file->Get("h_bco_diff_shifted");
    if (h_bco_diff_shifted == nullptr) {
      std::cerr << "ERROR: histogram h_bco_diff_shifted is missing in "
                << filename << std::endl;
      exit(1);
    }

    if (h_bco_diff_shifted->GetNbinsX() < observation.hist_bin) {
      std::cerr << "ERROR: histogram has fewer than " << observation.hist_bin
                << " bins in " << filename << std::endl;
      exit(1);
    }

    measured_n(row) = h_bco_diff_shifted->GetBinContent(observation.hist_bin);
    file->Close();
  }

  TVectorD solved_n;
  TMatrixD covariance_matrix;
  TMatrixD coefficient_matrix;

  const bool success = SolveShiftedSumWeighted(
      measured_n,
      measured_sigma,
      num_unknowns,
      equation_specs,
      offset_labels,
      &solved_n,
      &covariance_matrix,
      &coefficient_matrix);

  if (!success) {
    std::cerr << "Fit failed." << std::endl;
    return;
  }

  PrintSolvedEquations(num_unknowns, coefficient_matrix, measured_n, equation_specs, offset_labels);
  PrintSolution(num_unknowns, solved_n, covariance_matrix, offset_labels);
  SaveSolvedHistogram(num_unknowns, solved_n, covariance_matrix, offset_labels);
}






































int GetRunFromDelayValue(int delay) {
  std::vector< std::pair<int, int> > run_L1delay_list_scan6 = {
    {43291, 127},
    {43288, 126},
    {43285, 125},
    {43283, 124},
    {43282, 123},
    {43280, 122},
    {43278, 121},
    {43276, 120},
    {43313, 119}
  };
  std::vector< std::pair<int, int> > run_L1delay_list_scan7 = {
    {43408, 120},
    {43410, 119},
    {43412, 118},
    {43413, 117},
    {43414, 116},
    {43415, 115},
    {43417, 114},
    {43421, 113},
    {43426, 112},
    {43441, 110},
    {43434, 109},
    {43436, 108},
    {43438, 107},
    {43440, 106},
  };

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

int GetNumUnknowns(const std::vector<EquationSpec>& equation_specs) {
  int max_unknown_index = -1;
  for (const EquationSpec& equation : equation_specs) {
    for (int unknown_index : equation.unknown_indices) {
      max_unknown_index = std::max(max_unknown_index, unknown_index);
    }
  }

  if (max_unknown_index < 0) {
    std::cerr << "ERROR: No unknown terms were found in equation specs." << std::endl;
    exit(1);
  }

  return max_unknown_index + 1;
}

std::vector<EquationSpec> LoadEquationSpecs(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    std::cerr << "ERROR: Could not open equation spec file " << path << std::endl;
    exit(1);
  }

  const std::regex prefix_regex(
      R"(^N_\{relative_bin_position=(-?\d+), hist_bin=(\d+), delay=(\d+)\} = (.*)$)");
  const std::regex unknown_regex(R"(n_(\d+))");
  const std::regex offset_term_regex(
      R"(((?:\d+(?:\.\d+)?)\s*\*\s*)?(offset_[A-Za-z0-9_]+))");
  const std::regex offset_identity_regex(R"(offset_delay(\d+)_run(\d+).*)");

  std::vector<EquationSpec> equation_specs;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.find("No constraints here") != std::string::npos) {
      continue;
    }

    std::smatch prefix_match;
    if (!std::regex_match(line, prefix_match, prefix_regex)) {
      std::cerr << "ERROR: Could not parse equation line: " << line << std::endl;
      exit(1);
    }

    EquationSpec spec;
    spec.observation.relative_bin_position = std::stoi(prefix_match[1].str());
    spec.observation.hist_bin = std::stoi(prefix_match[2].str());
    spec.observation.delay = std::stoi(prefix_match[3].str());

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
      const std::string coefficient_text = (*it)[1].matched ? (*it)[1].str() : "";
      const std::string label = (*it)[2].str();

      double coefficient = 1.0;
      if (!coefficient_text.empty()) {
        std::string cleaned = coefficient_text;
        cleaned.erase(std::remove(cleaned.begin(), cleaned.end(), ' '), cleaned.end());
        if (!cleaned.empty() && cleaned.back() == '*') {
          cleaned.pop_back();
        }
        coefficient = std::stod(cleaned);
      }

      std::smatch offset_identity_match;
      if (std::regex_match(label, offset_identity_match, offset_identity_regex)) {
        const int offset_delay = std::stoi(offset_identity_match[1].str());
        const int offset_run = std::stoi(offset_identity_match[2].str());
        if (offset_delay != spec.observation.delay) {
          std::cerr << "ERROR: Delay mismatch in equation line: " << line << std::endl;
          exit(1);
        }
        if (GetRunFromDelayValue(offset_delay) != offset_run) {
          std::cerr << "ERROR: Run mismatch in equation line: " << line << std::endl;
          exit(1);
        }
      }

      spec.offset_terms.push_back({label, coefficient});
    }

    for (int unknown_index : spec.unknown_indices) {
      if (unknown_index < 0) {
        std::cerr << "ERROR: Unknown index out of range in equation line: " << line << std::endl;
        exit(1);
      }
    }

    // if (spec.observation.relative_bin_position < 1) {
    //   std::cerr << "ERROR: Unexpected relative_bin_position in equation line: " << line << std::endl;
    //   exit(1);
    // }
    if (spec.observation.hist_bin < 1) {
      std::cerr << "ERROR: Unexpected hist_bin in equation line: " << line << std::endl;
      exit(1);
    }

    equation_specs.push_back(spec);
  }

  if (equation_specs.empty()) {
    std::cerr << "ERROR: No constrained equations were loaded from " << path << std::endl;
    exit(1);
  }

  return equation_specs;
}

std::vector<std::string> CollectOffsetLabels(const std::vector<EquationSpec>& equation_specs) {
  std::vector<std::string> labels;
  std::map<std::string, int> label_to_index;

  for (const EquationSpec& equation : equation_specs) {
    for (const OffsetTerm& offset_term : equation.offset_terms) {
      if (label_to_index.find(offset_term.label) == label_to_index.end()) {
        label_to_index[offset_term.label] = static_cast<int>(labels.size());
        labels.push_back(offset_term.label);
      }
    }
  }

  return labels;
}

void BuildCoefficientMatrix(TMatrixD* matrix_a,
                            int num_unknowns,
                            const std::vector<EquationSpec>& equation_specs,
                            const std::vector<std::string>& offset_labels) {
  const int kNumParameters = num_unknowns + static_cast<int>(offset_labels.size());
  const int num_rows = static_cast<int>(equation_specs.size());

  std::map<std::string, int> label_to_index;
  for (int i = 0; i < static_cast<int>(offset_labels.size()); ++i) {
    label_to_index[offset_labels[i]] = i;
  }

  matrix_a->ResizeTo(num_rows, kNumParameters);
  matrix_a->Zero();

  for (int row = 0; row < num_rows; ++row) {
    const EquationSpec& equation = equation_specs[row];

    for (int unknown_index : equation.unknown_indices) {
      (*matrix_a)(row, unknown_index) = 1.0;
    }

    for (const OffsetTerm& offset_term : equation.offset_terms) {
      const auto it = label_to_index.find(offset_term.label);
      if (it == label_to_index.end()) {
        std::cerr << "ERROR: Unknown offset label in row " << row << std::endl;
        exit(1);
      }
      const int offset_col = num_unknowns + it->second;
      (*matrix_a)(row, offset_col) = offset_term.coefficient;
    }
  }
}

void SaveSolvedHistogram(int num_unknowns,
                         const TVectorD& solved_n,
                         const TMatrixD& covariance_matrix,
                         const std::vector<std::string>& offset_labels) {
  const char* output_dir =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
  const char* root_filename =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.root";
  const char* pdf_filename =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.pdf";

  gSystem->mkdir(output_dir, true);

  TH1D* h_solved_n = new TH1D("h_solved_n", "Solved n;bin index i;solved_n[i]",
                              num_unknowns, 0.5, num_unknowns + 0.5);
  for (int i = 0; i < num_unknowns; ++i) {
    h_solved_n->SetBinContent(i + 1, solved_n(i));
    if (covariance_matrix.GetNrows() > i && covariance_matrix.GetNcols() > i) {
      h_solved_n->SetBinError(i + 1, std::sqrt(covariance_matrix(i, i)));
    }
  }

  TH1D* h_solved_offsets = nullptr;
  if (!offset_labels.empty()) {
    h_solved_offsets =
        new TH1D("h_solved_offsets", "Solved offsets;offset label;value",
                 offset_labels.size(), 0.5, offset_labels.size() + 0.5);
    for (int i = 0; i < static_cast<int>(offset_labels.size()); ++i) {
      const int parameter_index = num_unknowns + i;
      h_solved_offsets->SetBinContent(i + 1, solved_n(parameter_index));
      if (covariance_matrix.GetNrows() > parameter_index &&
          covariance_matrix.GetNcols() > parameter_index) {
        h_solved_offsets->SetBinError(
            i + 1, std::sqrt(covariance_matrix(parameter_index, parameter_index)));
      }
      h_solved_offsets->GetXaxis()->SetBinLabel(i + 1, offset_labels[i].c_str());
    }
  }

  TFile output_file(root_filename, "RECREATE");
  h_solved_n->Write();
  if (h_solved_offsets != nullptr) {
    h_solved_offsets->Write();
  }
  output_file.Close();

  TCanvas canvas("c_solved_n", "c_solved_n", 900, 600);
  h_solved_n->SetStats(0);
  h_solved_n->Draw("HIST");
  canvas.SaveAs(pdf_filename);

  delete h_solved_n;
  delete h_solved_offsets;
}

bool SolveShiftedSumWeighted(const TVectorD& measured_n,
                             const TVectorD& measured_sigma,
                             int num_unknowns,
                             const std::vector<EquationSpec>& equation_specs,
                             const std::vector<std::string>& offset_labels,
                             TVectorD* solved_n,
                             TMatrixD* covariance_matrix,
                             TMatrixD* coefficient_matrix) {
  const int kNumParameters = num_unknowns + static_cast<int>(offset_labels.size());
  const int kNumObservations = static_cast<int>(equation_specs.size());

  if (measured_n.GetNrows() != kNumObservations ||
      measured_sigma.GetNrows() != kNumObservations) {
    std::cerr << "Input vector size mismatch." << std::endl;
    return false;
  }

  TMatrixD matrix_a;
  BuildCoefficientMatrix(&matrix_a, num_unknowns, equation_specs, offset_labels);
  if (coefficient_matrix != nullptr) {
    coefficient_matrix->ResizeTo(matrix_a);
    *coefficient_matrix = matrix_a;
  }

  TMatrixD weighted_matrix_a(kNumObservations, kNumParameters);
  TVectorD weighted_vector_y(kNumObservations);

  for (int row = 0; row < kNumObservations; ++row) {
    const double sigma = measured_sigma(row);
    if (sigma <= 0.0) {
      std::cerr << "Invalid sigma at row " << row << std::endl;
      return false;
    }

    weighted_vector_y(row) = measured_n(row) / sigma;
    for (int col = 0; col < kNumParameters; ++col) {
      weighted_matrix_a(row, col) = matrix_a(row, col) / sigma;
    }
  }

  TDecompSVD svd(weighted_matrix_a);
  Bool_t is_ok = kFALSE;

  solved_n->ResizeTo(kNumParameters);
  *solved_n = svd.Solve(weighted_vector_y, is_ok);

  if (!is_ok) {
    std::cerr << "SVD solve failed." << std::endl;
    return false;
  }

  TMatrixD matrix_a_transposed(TMatrixD::kTransposed, weighted_matrix_a);
  TMatrixD normal_matrix = matrix_a_transposed * weighted_matrix_a;

  TDecompSVD normal_svd(normal_matrix);
  Bool_t invert_ok = kFALSE;

  covariance_matrix->ResizeTo(kNumParameters, kNumParameters);
  *covariance_matrix = normal_svd.Invert(invert_ok);

  if (!invert_ok) {
    std::cerr << "Failed to invert normal matrix." << std::endl;
    return false;
  }

  return true;
}

void PrintSolvedEquations(int num_unknowns,
                          const TMatrixD& matrix_a,
                          const TVectorD& measured_n,
                          const std::vector<EquationSpec>& equation_specs,
                          const std::vector<std::string>& offset_labels) {
  const int kNumObservations = static_cast<int>(equation_specs.size());

  for (int row = 0; row < kNumObservations; ++row) {
    const ObservationIndex& observation = equation_specs[row].observation;
    std::cout << "N_{relative_bin_position=" << observation.relative_bin_position
              << ", hist_bin=" << observation.hist_bin
              << ", delay=" << observation.delay << "} = ";

    bool first_term = true;
    for (int col = 0; col < matrix_a.GetNcols(); ++col) {
      const double coefficient = matrix_a(row, col);
      if (std::abs(coefficient) < 1e-12) {
        continue;
      }

      if (!first_term) {
        std::cout << (coefficient >= 0.0 ? " + " : " - ");
      } else if (coefficient < 0.0) {
        std::cout << "-";
      }

      const double abs_coefficient = std::abs(coefficient);
      if (abs_coefficient != 1.0) {
        std::cout << abs_coefficient << " * ";
      }

      if (col < num_unknowns) {
        std::cout << "n_" << col + 1;
      } else {
        std::cout << offset_labels[col - num_unknowns];
      }
      first_term = false;
    }

    if (first_term) {
      std::cout << "0";
    }
    std::cout << " = " << measured_n(row) << std::endl;
  }
}

void PrintSolution(int num_unknowns,
                   const TVectorD& solved_n,
                   const TMatrixD& covariance_matrix,
                   const std::vector<std::string>& offset_labels) {
  for (int j = 1; j <= num_unknowns; ++j) {
    const double value = solved_n(j - 1);
    const double error = std::sqrt(covariance_matrix(j - 1, j - 1));
    std::cout << "n_" << j << " = " << value
              << " +/- " << error << std::endl;
  }

  for (int i = 0; i < static_cast<int>(offset_labels.size()); ++i) {
    const int parameter_index = num_unknowns + i;
    const double value = solved_n(parameter_index);
    const double error = std::sqrt(covariance_matrix(parameter_index, parameter_index));
    std::cout << offset_labels[i] << " = " << value
              << " +/- " << error << std::endl;
  }
}
