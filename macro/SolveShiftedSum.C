#include <cmath>
#include <fstream>
#include <iostream>
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
  int window_bin;
  int hist_bin;
  int delay;
};

struct EquationSpec {
  ObservationIndex observation;
  int offset_parity_index;
  std::vector<int> unknown_indices;
};

namespace {
constexpr int kDelayMin = 106;
constexpr int kDelayMax = 127;
constexpr int kUnknownBlockSize = 6;
constexpr int kNumUnknownBlocks = 6;
constexpr int kNumUnknowns = kUnknownBlockSize * kNumUnknownBlocks;
const char* kEquationSpecPath =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/macro/LinearEquationsToBeSolved.c";
}

int GetRunFromDelayValue(int delay);
int GetOffsetIndexFromDelay(int delay, int delay_min);
int GetDelayFromOffsetIndex(int offset_index, int delay_min);
int GetParityIndex(int bin);
const char* GetParityName(int parity_index);
int GetNumScannedPoints(int delay_min, int delay_max);
std::vector<EquationSpec> LoadEquationSpecs(const std::string& path);
void BuildCoefficientMatrix(TMatrixD* matrix_a, const std::vector<EquationSpec>& equation_specs, int delay_min, int delay_max);
void SaveSolvedHistogram(const TVectorD& solved_n);
bool SolveShiftedSumWeighted(const TVectorD& measured_n, const TVectorD& measured_sigma, const std::vector<EquationSpec>& equation_specs, TVectorD* solved_n, TMatrixD* covariance_matrix, TMatrixD* coefficient_matrix = nullptr);
void PrintSolvedEquations(const TMatrixD& matrix_a, const TVectorD& measured_n, const std::vector<EquationSpec>& equation_specs);
void PrintSolution(const TVectorD& solved_n, const TMatrixD& covariance_matrix);



void SolveShiftedSum() {
  const std::vector<EquationSpec> equation_specs = LoadEquationSpecs(kEquationSpecPath);
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
    std::string filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input/run%d.root", run);
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
      measured_n, measured_sigma, equation_specs,
      &solved_n, &covariance_matrix, &coefficient_matrix);

  if (!success) {
    std::cerr << "Fit failed." << std::endl;
    return;
  }

  PrintSolvedEquations(coefficient_matrix, measured_n, equation_specs);
  PrintSolution(solved_n, covariance_matrix);
  SaveSolvedHistogram(solved_n);
}


int GetRunFromDelayValue(int delay) {
  std::vector< pair<int, int> > run_L1delay_list_scan6 = {
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
  std::vector< pair<int, int> > run_L1delay_list_scan7 = {
    {43408, 120},
    {43410, 119},
    {43412, 118},
    {43413, 117},
    {43414, 116},
    {43415, 115},
    {43417, 114},
    {43421, 113},
    {43426, 112},
    // {43441, 111}, /* 111 does not exist ! later delete */
    {43441, 110},
    {43434, 109},
    {43436, 108},
    {43438, 107},
    {43440, 106},
  };

  for (int i=0; i<run_L1delay_list_scan6.size(); i++) {
    int run = run_L1delay_list_scan6[i].first;
    int l1_delay = run_L1delay_list_scan6[i].second;
    if (delay==l1_delay) {
      return run;
    }
  }

  for (int i=0; i<run_L1delay_list_scan7.size(); i++) {
    int run = run_L1delay_list_scan7[i].first;
    int l1_delay = run_L1delay_list_scan7[i].second;
    if (delay==l1_delay) {
      return run;
    }
  }
  return -1;
}

int GetNumScannedPoints(int delay_min, int delay_max) {
  return delay_max - delay_min + 1 - 1; /* -1 due to lack of a run with L1delay-111 */
}

int GetOffsetIndexFromDelay(int delay, int delay_min) {
  const int skipped_delay = delay_min + 5;
  if (delay < skipped_delay) {
    return delay - delay_min;
  }
  return delay - 1 - delay_min;
}

int GetDelayFromOffsetIndex(int offset_index, int delay_min) {
  const int skipped_delay = delay_min + 5;
  int delay = delay_min + offset_index;
  if (delay >= skipped_delay) {
    delay += 1;
  }
  return delay;
}

int GetParityIndex(int bin) {
  return (bin % 2 == 0) ? 1 : 0;
}

const char* GetParityName(int parity_index) {
  return (parity_index == 0) ? "odd" : "even";
}

std::vector<EquationSpec> LoadEquationSpecs(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    std::cerr << "ERROR: Could not open equation spec file " << path << std::endl;
    exit(1);
  }

  const std::regex prefix_regex(
      R"(^N_\{window_bin=(\d+), hist_bin=(\d+), delay=(\d+)\} = (.*)$)");
  const std::regex unknown_regex(R"(n_(\d+))");
  const std::regex offset_regex(R"(offset_delay(\d+)_run(\d+)_(odd|even))");

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
    spec.observation.window_bin = std::stoi(prefix_match[1].str());
    spec.observation.hist_bin = std::stoi(prefix_match[2].str());
    spec.observation.delay = std::stoi(prefix_match[3].str());

    const std::string rhs = prefix_match[4].str();
    for (std::sregex_iterator it(rhs.begin(), rhs.end(), unknown_regex), end; it != end; ++it) {
      spec.unknown_indices.push_back(std::stoi((*it)[1].str()) - 1);
    }

    std::smatch offset_match;
    if (!std::regex_search(rhs, offset_match, offset_regex)) {
      std::cerr << "ERROR: Missing offset term in equation line: " << line << std::endl;
      exit(1);
    }

    const int offset_delay = std::stoi(offset_match[1].str());
    const int offset_run = std::stoi(offset_match[2].str());
    const std::string parity_name = offset_match[3].str();
    spec.offset_parity_index = (parity_name == "odd") ? 0 : 1;

    if (offset_delay != spec.observation.delay) {
      std::cerr << "ERROR: Delay mismatch in equation line: " << line << std::endl;
      exit(1);
    }
    if (GetRunFromDelayValue(offset_delay) != offset_run) {
      std::cerr << "ERROR: Run mismatch in equation line: " << line << std::endl;
      exit(1);
    }
    if (spec.observation.window_bin < 1) {
      std::cerr << "ERROR: Unexpected window_bin in equation line: " << line << std::endl;
      exit(1);
    }
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

void BuildCoefficientMatrix(TMatrixD* matrix_a,
                            const std::vector<EquationSpec>& equation_specs,
                            int delay_min,
                            int delay_max) {
  const int num_scanned_points = GetNumScannedPoints(delay_min, delay_max);
  const int kNumOffsets = num_scanned_points * 2;
  const int kNumParameters = kNumUnknowns + kNumOffsets;
  const int num_rows = static_cast<int>(equation_specs.size());

  matrix_a->ResizeTo(num_rows, kNumParameters);
  matrix_a->Zero();

  for (int row = 0; row < num_rows; ++row) {
    const EquationSpec& equation = equation_specs[row];

    for (int unknown_index : equation.unknown_indices) {
      if (0 <= unknown_index && unknown_index < kNumUnknowns) {
        (*matrix_a)(row, unknown_index) = 1.0;
      } else {
        std::cerr << "ERROR: Unknown index out of range in row " << row << std::endl;
        exit(1);
      }
    }

    const int delay_offset_index = GetOffsetIndexFromDelay(equation.observation.delay, delay_min);
    const int offset_col = kNumUnknowns + 2 * delay_offset_index + equation.offset_parity_index;
    (*matrix_a)(row, offset_col) = 1.0;
  }
}

void SaveSolvedHistogram(const TVectorD& solved_n) {
  const char* output_dir =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
  const char* root_filename =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.root";
  const char* pdf_filename =
      "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output/solved_n.pdf";

  gSystem->mkdir(output_dir, true);

  TH1D* h_solved_n = new TH1D("h_solved_n", "Solved n;bin index i;solved_n[i]",
                              kNumUnknowns, 0.5, kNumUnknowns + 0.5);
  for (int i = 0; i < kNumUnknowns; ++i) {
    h_solved_n->SetBinContent(i + 1, solved_n(i));
  }

  TFile output_file(root_filename, "RECREATE");
  h_solved_n->Write();
  output_file.Close();

  TCanvas canvas("c_solved_n", "c_solved_n", 900, 600);
  h_solved_n->SetStats(0);
  h_solved_n->Draw("HIST");
  canvas.SaveAs(pdf_filename);

  delete h_solved_n;
}

bool SolveShiftedSumWeighted(const TVectorD& measured_n,
                             const TVectorD& measured_sigma,
                             const std::vector<EquationSpec>& equation_specs,
                             TVectorD* solved_n,
                             TMatrixD* covariance_matrix,
                             TMatrixD* coefficient_matrix = nullptr) {
  const int num_scanned_points = GetNumScannedPoints(kDelayMin, kDelayMax);
  const int kNumOffsets = num_scanned_points * 2;
  const int kNumParameters = kNumUnknowns + kNumOffsets;
  const int kNumObservations = static_cast<int>(equation_specs.size());

  if (measured_n.GetNrows() != kNumObservations ||
      measured_sigma.GetNrows() != kNumObservations) {
    std::cerr << "Input vector size mismatch." << std::endl;
    return false;
  }

  TMatrixD matrix_a;
  BuildCoefficientMatrix(&matrix_a, equation_specs, kDelayMin, kDelayMax);
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
  *covariance_matrix = normal_svd.Invert(invert_ok); /* (A^{T}*A)^{-1} */

  if (!invert_ok) {
    std::cerr << "Failed to invert normal matrix." << std::endl;
    return false;
  }

  /* Solution of the linear equations, that is (A^{T}*A)^{-1} * A^{T} * measured_n, will be calculated in  */
  return true;
}

void PrintSolvedEquations(const TMatrixD& matrix_a,
                          const TVectorD& measured_n,
                          const std::vector<EquationSpec>& equation_specs) {
  const int num_scanned_points = GetNumScannedPoints(kDelayMin, kDelayMax);
  const int kNumOffsets = num_scanned_points * 2;
  const int kNumObservations = static_cast<int>(equation_specs.size());

  for (int row = 0; row < kNumObservations; ++row) {
    const ObservationIndex& observation = equation_specs[row].observation;
    std::cout << "N_{window_bin=" << observation.window_bin
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

      if (col < kNumUnknowns) {
        std::cout << "n_" << col + 1;
      } else if (col < kNumUnknowns + kNumOffsets) {
        const int local_offset_index = col - kNumUnknowns;
        const int offset_delay = GetDelayFromOffsetIndex(local_offset_index / 2, kDelayMin);
        const int parity_index = local_offset_index % 2;
        const int run = GetRunFromDelayValue(offset_delay);
        std::cout << "offset_delay" << offset_delay << "_run" << run
                  << "_" << GetParityName(parity_index);
      } else {
        std::cout << "parameter_" << col + 1;
      }
      first_term = false;
    }

    if (first_term) {
      std::cout << "0";
    }
    std::cout << " = " << measured_n(row) << std::endl;
  }
}

void PrintSolution(const TVectorD& solved_n, const TMatrixD& covariance_matrix) {
  const int num_scanned_points = GetNumScannedPoints(kDelayMin, kDelayMax);
  const int kNumOffsets = num_scanned_points * 2;

  /* Fine timing distribution part */
  for (int j = 1; j <= kNumUnknowns; ++j) {
    const double value = solved_n(j - 1);
    const double error = std::sqrt(covariance_matrix(j - 1, j - 1));
    std::cout << "n_" << j << " = " << value
              << " +/- " << error << std::endl;
  }

  /* Offset part */
  for (int i = 0; i < kNumOffsets; ++i) {
    const int delay = GetDelayFromOffsetIndex(i / 2, kDelayMin);
    const int parity_index = i % 2;
    const int run = GetRunFromDelayValue(delay);
    const double offset = solved_n(kNumUnknowns + i);
    const double offset_error = std::sqrt(covariance_matrix(kNumUnknowns + i, kNumUnknowns + i));
    std::cout << "offset_delay" << delay << "_run" << run
              << "_" << GetParityName(parity_index)
              << " = " << offset
              << " +/- " << offset_error << std::endl;
  }
}
