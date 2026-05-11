#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"

namespace {

const char* kDefaultInputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input";
const char* kDefaultOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input";
const char* kTreeName = "tree";
const char* kCutBranchName = "n_hits_around_trigger_timing";

std::string BuildInputPath(int run) {
  return Form("%s/run%d.root", kDefaultInputDir, run);
}

uint16_t CalculateBcoDiff(unsigned long bco_full, unsigned short fphx_bco) {
  const uint16_t bco_low7 = static_cast<uint16_t>(bco_full & 0x7fUL);
  const uint16_t fphx_bco_u16 = static_cast<uint16_t>(fphx_bco);

  if (bco_low7 > fphx_bco_u16) {
    return static_cast<uint16_t>((fphx_bco_u16 + 128U) - bco_low7);
  }
  return static_cast<uint16_t>(fphx_bco_u16 - bco_low7);
}

}  // namespace

void MakeMultiplicityCutBcoDiffHistogram(
    int run = 43276,
    int multiplicity_min = 10,
    int multiplicity_max = 100,
    bool normalize_per_selected_event = true) {
  if (multiplicity_min >= multiplicity_max) {
    std::cerr << "ERROR: multiplicity_min must be smaller than multiplicity_max." << std::endl;
    return;
  }

  const std::string input_path = BuildInputPath(run);
  std::unique_ptr<TFile> input_file(TFile::Open(input_path.c_str(), "READ"));
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open input file: " << input_path << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(input_file->Get(kTreeName));
  if (tree == nullptr) {
    std::cerr << "ERROR: " << kTreeName << " is missing in " << input_path << std::endl;
    return;
  }

  if (tree->GetBranch("bco") == nullptr ||
      tree->GetBranch("FPHX_BCO") == nullptr ||
      tree->GetBranch(kCutBranchName) == nullptr) {
    std::cerr << "ERROR: bco, FPHX_BCO, and/or " << kCutBranchName
              << " branch is missing in " << input_path << std::endl;
    return;
  }

  std::vector<unsigned long>* bco = nullptr;
  std::vector<unsigned short>* fphx_bco = nullptr;
  int cut_value = 0;

  tree->SetBranchAddress("bco", &bco);
  tree->SetBranchAddress("FPHX_BCO", &fphx_bco);
  tree->SetBranchAddress(kCutBranchName, &cut_value);

  const TString hist_name = Form("h_bco_diff_%d_lt_%s_lt_%d",
                                 multiplicity_min,
                                 kCutBranchName,
                                 multiplicity_max);
  TH1D hist(hist_name,
            Form("%s;BCO diff;INTT hit counts per selected event", hist_name.Data()),
            128,
            0.0,
            128.0);
  hist.Sumw2();

  Long64_t selected_events = 0;
  Long64_t truncated_entries = 0;
  const Long64_t num_entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < num_entries; ++entry) {
    tree->GetEntry(entry);

    if (!(multiplicity_min < cut_value && cut_value < multiplicity_max)) {
      continue;
    }

    ++selected_events;

    if (bco == nullptr || fphx_bco == nullptr) {
      std::cerr << "WARNING: Null vector branch pointer at entry " << entry << std::endl;
      continue;
    }

    const std::size_t num_hits = std::min(bco->size(), fphx_bco->size());
    if (bco->size() != fphx_bco->size()) {
      ++truncated_entries;
    }

    for (std::size_t index = 0; index < num_hits; ++index) {
      hist.Fill(CalculateBcoDiff(bco->at(index), fphx_bco->at(index)));
    }
  }

  if (normalize_per_selected_event && selected_events > 0) {
    hist.Scale(1.0 / static_cast<double>(selected_events));
  }

  const TString output_path = Form("%s/run%d_multiplicity_cut.root", kDefaultOutputDir, run);
  std::unique_ptr<TFile> output_file(TFile::Open(output_path, "RECREATE"));
  if (!output_file || output_file->IsZombie()) {
    std::cerr << "ERROR: Could not create output file: " << output_path << std::endl;
    return;
  }

  output_file->cd();
  hist.Write();
  output_file->Close();

  std::cout << "Created " << output_path << std::endl;
  std::cout << "Histogram name: " << hist_name << std::endl;
  std::cout << "Applied cut: " << multiplicity_min << " < " << kCutBranchName
            << " < " << multiplicity_max << std::endl;
  std::cout << "Selected events: " << selected_events << " / " << num_entries << std::endl;
  if (truncated_entries > 0) {
    std::cout << "WARNING: Truncated " << truncated_entries
              << " entries due to vector size mismatches between bco and FPHX_BCO." << std::endl;
  }
}
