#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
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
const unsigned int kFirstFelixServer = 3001;
const int kNumFelixServers = 8;
const int kNumFelixChannels = 14;
const int kNumHalfLadders = kNumFelixServers * kNumFelixChannels;

std::string BuildInputPath(int run) {
  return Form("%s/run%d.root", kDefaultInputDir, run);
}

std::string BuildCutSuffix(int multiplicity_min, int multiplicity_max) {
  std::ostringstream suffix;
  suffix << multiplicity_min << "to" << multiplicity_max;
  return suffix.str();
}

TH1D* MakeHistogram(const std::string& name) {
  TH1D* hist = new TH1D(name.c_str(),
                        Form("%s;BCO diff;INTT hit counts per selected event", name.c_str()),
                        128,
                        0.0,
                        128.0);
  hist->Sumw2();
  return hist;
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
      tree->GetBranch("felix_server") == nullptr ||
      tree->GetBranch("felix_channel") == nullptr ||
      tree->GetBranch(kCutBranchName) == nullptr) {
    std::cerr << "ERROR: bco, FPHX_BCO, felix_server, felix_channel, and/or " << kCutBranchName
              << " branch is missing in " << input_path << std::endl;
    return;
  }

  std::vector<unsigned long>* bco = nullptr;
  std::vector<unsigned short>* fphx_bco = nullptr;
  std::vector<unsigned int>* felix_server = nullptr;
  std::vector<unsigned short>* felix_channel = nullptr;
  int cut_value = 0;

  tree->SetBranchAddress("bco", &bco);
  tree->SetBranchAddress("FPHX_BCO", &fphx_bco);
  tree->SetBranchAddress("felix_server", &felix_server);
  tree->SetBranchAddress("felix_channel", &felix_channel);
  tree->SetBranchAddress(kCutBranchName, &cut_value);

  const std::string cut_suffix = BuildCutSuffix(multiplicity_min, multiplicity_max);
  const std::string hist_name = "h_bco_diff_" + cut_suffix;
  std::unique_ptr<TH1D> hist(MakeHistogram(hist_name));

  std::vector<std::unique_ptr<TH1D>> felix_server_hists;
  felix_server_hists.reserve(kNumFelixServers);
  for (int server_index = 0; server_index < kNumFelixServers; ++server_index) {
    felix_server_hists.emplace_back(
        MakeHistogram(Form("h_bco_diff_%s_felix_server_%d", cut_suffix.c_str(), server_index)));
  }

  std::vector<std::unique_ptr<TH1D>> halfladder_hists;
  halfladder_hists.reserve(kNumHalfLadders);
  for (int halfladder_id = 0; halfladder_id < kNumHalfLadders; ++halfladder_id) {
    halfladder_hists.emplace_back(
        MakeHistogram(Form("h_bco_diff_%s_halfladder_id_%d", cut_suffix.c_str(), halfladder_id)));
  }

  Long64_t selected_events = 0;
  Long64_t truncated_entries = 0;
  Long64_t skipped_hits_due_to_invalid_mapping = 0;
  const Long64_t num_entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < num_entries; ++entry) {
    tree->GetEntry(entry);

    if (!(multiplicity_min < cut_value && cut_value < multiplicity_max)) {
      continue;
    }

    ++selected_events;

    if (bco == nullptr || fphx_bco == nullptr || felix_server == nullptr || felix_channel == nullptr) {
      std::cerr << "WARNING: Null vector branch pointer at entry " << entry << std::endl;
      continue;
    }

    const std::size_t num_hits =
        std::min(std::min(bco->size(), fphx_bco->size()),
                 std::min(felix_server->size(), felix_channel->size()));
    if (bco->size() != fphx_bco->size() ||
        bco->size() != felix_server->size() ||
        bco->size() != felix_channel->size()) {
      ++truncated_entries;
    }

    for (std::size_t index = 0; index < num_hits; ++index) {
      const uint16_t bco_diff = CalculateBcoDiff(bco->at(index), fphx_bco->at(index));
      hist->Fill(bco_diff);

      const int server_index = static_cast<int>(felix_server->at(index)) - static_cast<int>(kFirstFelixServer);
      if (server_index < 0 || server_index >= kNumFelixServers) {
        ++skipped_hits_due_to_invalid_mapping;
        continue;
      }
      felix_server_hists[server_index]->Fill(bco_diff);

      const int channel = static_cast<int>(felix_channel->at(index));
      if (channel < 0 || channel >= kNumFelixChannels) {
        ++skipped_hits_due_to_invalid_mapping;
        continue;
      }

      const int halfladder_id = kNumFelixChannels * server_index + channel;
      halfladder_hists[halfladder_id]->Fill(bco_diff);
    }
  }

  if (normalize_per_selected_event && selected_events > 0) {
    const double scale = 1.0 / static_cast<double>(selected_events);
    hist->Scale(scale);
    for (std::size_t server_index = 0; server_index < felix_server_hists.size(); ++server_index) {
      felix_server_hists[server_index]->Scale(scale);
    }
    for (std::size_t halfladder_id = 0; halfladder_id < halfladder_hists.size(); ++halfladder_id) {
      halfladder_hists[halfladder_id]->Scale(scale);
    }
  }

  const TString output_path = Form("%s/run%d_multiplicity_%dto%d.root", kDefaultOutputDir, run, multiplicity_min, multiplicity_max);
  std::unique_ptr<TFile> output_file(TFile::Open(output_path, "RECREATE"));
  if (!output_file || output_file->IsZombie()) {
    std::cerr << "ERROR: Could not create output file: " << output_path << std::endl;
    return;
  }

  output_file->cd();
  hist->Write();
  for (std::size_t server_index = 0; server_index < felix_server_hists.size(); ++server_index) {
    felix_server_hists[server_index]->Write();
  }
  for (std::size_t halfladder_id = 0; halfladder_id < halfladder_hists.size(); ++halfladder_id) {
    halfladder_hists[halfladder_id]->Write();
  }
  output_file->Close();

  std::cout << "Created " << output_path << std::endl;
  std::cout << "Histogram name: " << hist_name << std::endl;
  std::cout << "Applied cut: " << multiplicity_min << " < " << kCutBranchName
            << " < " << multiplicity_max << std::endl;
  std::cout << "Selected events: " << selected_events << " / " << num_entries << std::endl;
  if (truncated_entries > 0) {
    std::cout << "WARNING: Truncated " << truncated_entries
              << " entries due to vector size mismatches between hit branches." << std::endl;
  }
  if (skipped_hits_due_to_invalid_mapping > 0) {
    std::cout << "WARNING: Skipped " << skipped_hits_due_to_invalid_mapping
              << " hits due to invalid felix_server/felix_channel values." << std::endl;
  }
}
