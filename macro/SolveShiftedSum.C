#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TSystem.h"
#include "TVectorD.h"

struct ObservationIndex {
  int i;
  int delay; /* [1/6 BCO] */
};

int GetRunFromDelayValue(int delay);
int GetOffsetIndexFromDelay(int delay, int delay_min);
int GetDelayFromOffsetIndex(int offset_index, int delay_min);
int GetParityIndex(int bin);
const char* GetParityName(int parity_index);

int GetRowIndex(int i, int delay, int delay_min, int num_scanned_points) {
  if (delay < 5){ /* L1delay 111 = 106+5*/
    return (i - 1) * num_scanned_points + (delay - delay_min);
  } else {
    /* Since a run with L1delay=111 is absent */
    return (i - 1) * num_scanned_points + (delay-1 - delay_min);
  }
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

void BuildCoefficientMatrix(TMatrixD* matrix_a, int d_min, int d_max) {
  const int kNumUnknowns = 42;
  const int kNumGroups = 7;
  const int num_scanned_points = d_max - d_min + 1 -1; /* -1 due to lack of a run with L1delay-111 */
  const int kNumOffsets = num_scanned_points * 2;
  const int kNumParameters = kNumUnknowns + kNumOffsets;
  const int num_rows = kNumGroups * num_scanned_points;

  matrix_a->ResizeTo(num_rows, kNumParameters);
  matrix_a->Zero();

  for (int d = d_min; d <= d_max; ++d) {
    if (d==5) {continue;}
    for (int i = 1; i <= kNumGroups; ++i) {
      const int row = GetRowIndex(i, d, d_min, num_scanned_points);

      const int j_start = 6 * i + d - 5;
      const int j_end = 6 * i + d;

      for (int j = j_start; j <= j_end; ++j) {
        if (1 <= j && j <= kNumUnknowns) {
          (*matrix_a)(row, j - 1) = 1.0;
        }
      }

      const int delay_offset_index = GetOffsetIndexFromDelay(d, d_min);
      const int parity_index = GetParityIndex(i);
      const int offset_col = kNumUnknowns + 2 * delay_offset_index + parity_index;
      (*matrix_a)(row, offset_col) = 1.0;
    }
  }
}

void SaveSolvedHistogram(const TVectorD& solved_n) {
  const int kNumUnknowns = 42;
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
                             TVectorD* solved_n,
                             TMatrixD* covariance_matrix,
                             TMatrixD* coefficient_matrix = nullptr) {
  const int kNumUnknowns = 42; /* 7 times 6 */
  const int kDelayMin = 106;
  const int kDelayMax = 127;
  const int kNumGroups = 7;
  const int num_scanned_points = kDelayMax - kDelayMin + 1 -1; /* -1 due to lack of a run with L1delay-111 */
  const int kNumOffsets = num_scanned_points * 2;
  const int kNumParameters = kNumUnknowns + kNumOffsets;
  const int kNumObservations = kNumGroups * num_scanned_points;

  if (measured_n.GetNrows() != kNumObservations ||
      measured_sigma.GetNrows() != kNumObservations) {
    std::cerr << "Input vector size mismatch." << std::endl;
    return false;
  }

  TMatrixD matrix_a;
  BuildCoefficientMatrix(&matrix_a, kDelayMin-106, kDelayMax-106);
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

void PrintSolvedEquations(const TMatrixD& matrix_a, const TVectorD& measured_n) {
  const int kNumUnknowns = 42;
  const int kDelayMin = 106;
  const int kDelayMax = 127;
  const int kNumGroups = 7;
  const int kNumOffsets = (kDelayMax - kDelayMin + 1 - 1) * 2;
  const int num_scanned_points = kDelayMax - kDelayMin + 1 - 1; /* -1 due to lack of a run with L1delay-111 */

  for (int delay = kDelayMin; delay <= kDelayMax; ++delay) {
    if (delay == 111) {
      continue;
    }

    const int d = delay - kDelayMin;
    for (int i = 1; i <= kNumGroups; ++i) {
      const int row = GetRowIndex(i, d, 0, num_scanned_points);
      std::cout << "N_" << i << "^{delay=" << delay << "} = ";

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
}

void PrintSolution(const TVectorD& solved_n, const TMatrixD& covariance_matrix) {
  const int kNumUnknowns = 42;
  const int kDelayMin = 106;
  const int kDelayMax = 127;
  const int kNumOffsets = (kDelayMax - kDelayMin + 1 - 1) * 2;

  for (int j = 1; j <= kNumUnknowns; ++j) {
    const double value = solved_n(j - 1);
    const double error = std::sqrt(covariance_matrix(j - 1, j - 1));
    std::cout << "n_" << j << " = " << value
              << " +/- " << error << std::endl;
  }

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

void SolveShiftedSum() {
  const int kDelayMin = 106;
  const int kDelayMax = 127;
  const int kNumGroups = 7;
  const int num_scanned_points = kDelayMax - kDelayMin + 1 -1; /* -1 due to lack of a run with L1delay-111 */
  const int kNumObservations = kNumGroups * num_scanned_points;

  TVectorD measured_n(kNumObservations);
  TVectorD measured_sigma(kNumObservations);

  measured_n.Zero();
  measured_sigma.Zero();

  // 例: 全点の誤差を 1.0 にしておく
  for (int row = 0; row < kNumObservations; ++row) {
    measured_sigma(row) = 0.01;
  }

  // 実データ代入
  for (int delay=kDelayMin; delay<=kDelayMax; delay++) {
    if (delay==111) {
      continue; /* skip L1delay=111 since we don't have such run */
    }
    int run = GetRunFromDelayValue(delay);
    std::string filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input/run%d.root", run);
    TFile* file = TFile::Open(filename.c_str(), "update");
    if ( !file || file->IsZombie() ){
      std::cerr << "ERROR: Could not open" << filename << std::endl;
      exit(1);
    }
    for (int bin=1; bin <= 7; bin++) {
      // TH1D* h_bco_diff_shifted = (TH1D*)file->Get("h_bco_diff_shifted");
      TH1D* h_bco_diff_shifted = (TH1D*)file->Get("h_bco_diff_shifted_minus_plateau");
      measured_n(GetRowIndex(bin, delay-106, kDelayMin-106, num_scanned_points)) = h_bco_diff_shifted->GetBinContent(bin);
    }
  }

  TVectorD solved_n;
  TMatrixD covariance_matrix;
  TMatrixD coefficient_matrix;

  const bool success = SolveShiftedSumWeighted(measured_n, measured_sigma, &solved_n, &covariance_matrix, &coefficient_matrix);

  if (!success) {
    std::cerr << "Fit failed." << std::endl;
    return;
  }

  PrintSolvedEquations(coefficient_matrix, measured_n);
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
