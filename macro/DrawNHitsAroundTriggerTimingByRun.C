#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <dirent.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"

namespace {

const char* kDefaultInputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input";
const char* kDefaultOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";
const char* kTreeName = "tree";
const char* kBranchName = "n_hits_around_trigger_timing";
const int kXAxisMin = 0;
const int kXAxisMax = 1000;

struct RunFileInfo {
  int run = -1;
  std::string filepath;
};

struct RunGroup {
  std::string label;
  std::vector<int> runs;
};

struct RunHistogram {
  int run = -1;
  TH1D* hist = nullptr;
};

bool IsTargetFile(const char* name) {
  const TString filename(name);
  return filename.BeginsWith("run") && filename.EndsWith(".root");
}

int ExtractRunNumber(const std::string& filepath) {
  const std::size_t run_pos = filepath.rfind("run");
  const std::size_t dot_pos = filepath.rfind(".root");
  if (run_pos == std::string::npos || dot_pos == std::string::npos || dot_pos <= run_pos + 3) {
    return -1;
  }
  return std::stoi(filepath.substr(run_pos + 3, dot_pos - (run_pos + 3)));
}

std::vector<int> BuildLineColors() {
  return {
      kRed + 1,
      kBlue + 1,
      kGreen + 2,
      kMagenta + 1,
      kOrange + 7,
      kCyan + 2,
      kPink + 7,
      kSpring + 5,
      kViolet + 2,
      kAzure + 2,
      kTeal + 3,
      kYellow + 2,
  };
}

std::vector<RunGroup> BuildRunGroups(const std::vector<RunFileInfo>& run_files) {
  std::vector<RunGroup> groups = {
      {"Group1", {43291, 43288, 43285, 43283, 43282, 43280, 43278, 43276, 43313}},
      {"Group2", {43408, 43410, 43412, 43413, 43414, 43415, 43417, 43421, 43426}},
  };

  std::set<int> assigned_runs;
  for (const RunGroup& group : groups) {
    assigned_runs.insert(group.runs.begin(), group.runs.end());
  }

  RunGroup group3;
  group3.label = "Group3";
  for (const RunFileInfo& run_file : run_files) {
    if (assigned_runs.count(run_file.run) == 0) {
      group3.runs.push_back(run_file.run);
    }
  }
  groups.push_back(group3);

  return groups;
}

const RunHistogram* FindRunHistogram(const std::vector<RunHistogram>& run_histograms, int run) {
  for (const RunHistogram& run_histogram : run_histograms) {
    if (run_histogram.run == run) {
      return &run_histogram;
    }
  }
  return nullptr;
}

bool CollectRunFiles(const std::string& input_dir, std::vector<RunFileInfo>* run_files) {
  if (run_files == nullptr) {
    return false;
  }

  DIR* dir = opendir(input_dir.c_str());
  if (!dir) {
    std::cerr << "ERROR: Could not open directory " << input_dir << std::endl;
    return false;
  }

  std::vector<std::string> filepaths;
  while (dirent* entry = readdir(dir)) {
    if (!IsTargetFile(entry->d_name)) {
      continue;
    }
    filepaths.push_back(input_dir + "/" + entry->d_name);
  }
  closedir(dir);

  std::sort(filepaths.begin(), filepaths.end());
  if (filepaths.empty()) {
    std::cerr << "ERROR: No run*.root files found in " << input_dir << std::endl;
    return false;
  }

  run_files->clear();
  run_files->reserve(filepaths.size());

  for (const std::string& filepath : filepaths) {
    const int run = ExtractRunNumber(filepath);
    if (run < 0) {
      std::cerr << "WARNING: Could not parse run number from " << filepath << std::endl;
      continue;
    }

    std::unique_ptr<TFile> input_file(TFile::Open(filepath.c_str(), "READ"));
    if (!input_file || input_file->IsZombie()) {
      std::cerr << "WARNING: Could not open " << filepath << std::endl;
      continue;
    }

    TTree* tree = dynamic_cast<TTree*>(input_file->Get(kTreeName));
    if (tree == nullptr) {
      std::cerr << "WARNING: " << kTreeName << " is missing in " << filepath << std::endl;
      continue;
    }
    if (tree->GetBranch(kBranchName) == nullptr) {
      std::cerr << "WARNING: " << kBranchName << " is missing in " << filepath << std::endl;
      continue;
    }

    RunFileInfo info;
    info.run = run;
    info.filepath = filepath;
    run_files->push_back(info);
  }

  if (run_files->empty()) {
    std::cerr << "ERROR: No usable ROOT files found in " << input_dir << std::endl;
    return false;
  }

  std::sort(run_files->begin(), run_files->end(),
            [](const RunFileInfo& left, const RunFileInfo& right) { return left.run < right.run; });
  return true;
}

TH1D* BuildHistogram(const RunFileInfo& run_file, bool normalize) {
  std::unique_ptr<TFile> input_file(TFile::Open(run_file.filepath.c_str(), "READ"));
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << run_file.filepath << std::endl;
    return nullptr;
  }

  TTree* tree = dynamic_cast<TTree*>(input_file->Get(kTreeName));
  if (tree == nullptr) {
    std::cerr << "ERROR: " << kTreeName << " is missing in " << run_file.filepath << std::endl;
    return nullptr;
  }

  int n_hits_around_trigger_timing = 0;
  tree->SetBranchAddress(kBranchName, &n_hits_around_trigger_timing);

  TH1D* hist = new TH1D(Form("h_%s_run%d", kBranchName, run_file.run),
                        "",
                        kXAxisMax - kXAxisMin + 1,
                        kXAxisMin - 0.5,
                        kXAxisMax + 0.5);
  hist->SetDirectory(nullptr);

  const Long64_t num_entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < num_entries; ++entry) {
    tree->GetEntry(entry);
    hist->Fill(n_hits_around_trigger_timing);
  }

  if (normalize && hist->Integral() > 0.0) {
    hist->Scale(1.0 / hist->Integral());
  }

  return hist;
}

bool BuildRunHistograms(const std::vector<RunFileInfo>& run_files,
                        bool normalize,
                        std::vector<RunHistogram>* run_histograms,
                        double* common_y_max) {
  if (run_histograms == nullptr || common_y_max == nullptr) {
    return false;
  }

  run_histograms->clear();
  run_histograms->reserve(run_files.size());
  *common_y_max = 0.0;

  for (const RunFileInfo& run_file : run_files) {
    TH1D* hist = BuildHistogram(run_file, normalize);
    if (hist == nullptr) {
      for (const RunHistogram& run_histogram : *run_histograms) {
        delete run_histogram.hist;
      }
      run_histograms->clear();
      return false;
    }

    RunHistogram run_histogram;
    run_histogram.run = run_file.run;
    run_histogram.hist = hist;
    run_histograms->push_back(run_histogram);
    *common_y_max = std::max(*common_y_max, hist->GetMaximum());
  }

  return !run_histograms->empty();
}

bool DrawRunComparison(const std::string& canvas_name,
                       const std::string& canvas_title,
                       const std::vector<int>& runs,
                       const std::vector<RunHistogram>& run_histograms,
                       const std::string& output_path,
                       const std::vector<int>& line_colors,
                       bool normalize,
                       double common_y_max) {
  std::vector<const RunHistogram*> selected_histograms;
  selected_histograms.reserve(runs.size());

  for (int run : runs) {
    const RunHistogram* run_histogram = FindRunHistogram(run_histograms, run);
    if (run_histogram == nullptr) {
      std::cerr << "WARNING: Run " << run << " is not available and will be skipped." << std::endl;
      continue;
    }
    selected_histograms.push_back(run_histogram);
  }

  if (selected_histograms.empty()) {
    std::cerr << "WARNING: No histograms were selected for " << canvas_name << std::endl;
    return false;
  }

  TCanvas* canvas = new TCanvas(canvas_name.c_str(), canvas_title.c_str(), 1400, 900);
  canvas->SetMargin(0.10, 0.04, 0.12, 0.08);
  canvas->SetTicks(1, 1);
  canvas->SetLogy(true);

  TLegend* legend = new TLegend(0.12, 0.74, 0.88, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetNColumns(std::min(static_cast<int>(selected_histograms.size()), 4));
  legend->SetTextSize(0.028);

  const double y_min = normalize ? 1e-6 : 0.5;
  const double y_max = std::max(common_y_max * 2.0, y_min * 10.0);

  for (std::size_t index = 0; index < selected_histograms.size(); ++index) {
    TH1D* hist = selected_histograms[index]->hist;
    hist->SetLineColor(line_colors[index % line_colors.size()]);
    hist->SetLineWidth(2);
    hist->SetTitle(canvas_title.c_str());
    hist->GetXaxis()->SetTitle("n_hits_around_trigger_timing");
    hist->GetXaxis()->SetRangeUser(kXAxisMin, kXAxisMax);
    hist->GetYaxis()->SetTitle(normalize ? "Normalized entries" : "Entries");
    hist->GetYaxis()->SetRangeUser(y_min, y_max);

    if (index == 0) {
      hist->Draw("hist");
    } else {
      hist->Draw("hist same");
    }

    legend->AddEntry(hist, Form("Run %d", selected_histograms[index]->run), "l");
  }

  legend->Draw();
  canvas->Print(output_path.c_str());
  std::cout << "Saved " << output_path << std::endl;

  delete legend;
  delete canvas;
  return true;
}

}  // namespace

void DrawNHitsAroundTriggerTimingByRun(
    const char* input_dir = kDefaultInputDir,
    const char* output_dir = kDefaultOutputDir,
    bool normalize = true) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(output_dir, true);

  std::vector<RunFileInfo> run_files;
  if (!CollectRunFiles(input_dir, &run_files)) {
    return;
  }

  const std::vector<int> line_colors = BuildLineColors();
  std::vector<RunHistogram> run_histograms;
  double common_y_max = 0.0;
  if (!BuildRunHistograms(run_files, normalize, &run_histograms, &common_y_max)) {
    return;
  }

  std::vector<int> all_runs;
  all_runs.reserve(run_histograms.size());
  for (const RunHistogram& run_histogram : run_histograms) {
    all_runs.push_back(run_histogram.run);
  }

  const TString suffix = normalize ? "normalized" : "raw";
  int num_success = 0;
  if (DrawRunComparison("c_n_hits_around_trigger_timing_by_run_all",
                        "n_hits_around_trigger_timing comparison by run",
                        all_runs,
                        run_histograms,
                        Form("%s/n_hits_around_trigger_timing_by_run_%s.pdf",
                             output_dir,
                             suffix.Data()),
                        line_colors,
                        normalize,
                        common_y_max)) {
    ++num_success;
  }

  const std::vector<RunGroup> groups = BuildRunGroups(run_files);
  for (const RunGroup& group : groups) {
    if (group.runs.empty()) {
      std::cerr << "WARNING: " << group.label << " has no runs to draw." << std::endl;
      continue;
    }

    if (DrawRunComparison(Form("c_n_hits_around_trigger_timing_by_run_%s", group.label.c_str()),
                          Form("n_hits_around_trigger_timing %s", group.label.c_str()),
                          group.runs,
                          run_histograms,
                          Form("%s/n_hits_around_trigger_timing_by_run_%s_%s.pdf",
                               output_dir,
                               group.label.c_str(),
                               suffix.Data()),
                          line_colors,
                          normalize,
                          common_y_max)) {
      ++num_success;
    }
  }

  std::cout << "Created " << num_success << " comparison canvases." << std::endl;

  for (const RunHistogram& run_histogram : run_histograms) {
    delete run_histogram.hist;
  }
}
