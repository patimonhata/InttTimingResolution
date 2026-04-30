#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

namespace {

struct RunGroup {
  const char* label;
  std::vector<int> runs;
};

const char* kDefaultInputDir =
    // "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input/AsOf2024";
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/input";
const char* kDefaultOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";

std::vector<RunGroup> BuildRunGroups() {
  return {
      {"Group1", {43291, 43288, 43285, 43283, 43282, 43280, 43278, 43276, 43313}},
      {"Group2", {43408, 43410, 43412, 43413, 43414, 43415, 43417, 43421, 43426}},
      {"Group3", {43441, 43436, 43438, 43440}},
  };
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
  };
}

TH1D* LoadShiftedHistogram(const std::string& input_dir, int run) {
  const TString filepath = Form("%s/run%d.root", input_dir.c_str(), run);
  TFile* input_file = TFile::Open(filepath, "READ");
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << filepath << std::endl;
    delete input_file;
    return nullptr;
  }

  TH1D* input_hist = dynamic_cast<TH1D*>(input_file->Get("h_bco_diff_shifted"));
  if (input_hist == nullptr) {
    std::cerr << "ERROR: h_bco_diff_shifted is missing in " << filepath << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(input_hist->Clone(Form("h_bco_diff_shifted_run%d", run)));
  if (hist == nullptr) {
    std::cerr << "ERROR: Failed to clone h_bco_diff_shifted from " << filepath << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  hist->SetDirectory(nullptr);
  input_file->Close();
  delete input_file;
  return hist;
}

bool DrawGroup(const RunGroup& group,
               const std::string& input_dir,
               const std::string& output_dir,
               const std::vector<int>& line_colors) {
  TCanvas* canvas = new TCanvas(Form("c_%s", group.label), group.label, 900, 700);
  canvas->SetMargin(0.12, 0.04, 0.11, 0.06);
  canvas->SetTicks(1, 1);
  canvas->SetLogy(true);

  TLegend* legend = new TLegend(0.62, 0.60, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);

  std::vector<TH1D*> histograms;
  histograms.reserve(group.runs.size());

  bool first_hist = true;
  for (std::size_t index = 0; index < group.runs.size(); ++index) {
    const int run = group.runs[index];
    TH1D* hist = LoadShiftedHistogram(input_dir, run);
    if (hist == nullptr) {
      for (TH1D* loaded_hist : histograms) {
        delete loaded_hist;
      }
      delete legend;
      delete canvas;
      return false;
    }

    hist->SetLineColor(line_colors[index % line_colors.size()]);
    hist->SetLineWidth(2);
    hist->SetStats(0);
    hist->SetTitle(Form("%s h_bco_diff_shifted comparison", group.label));
    hist->GetYaxis()->SetRangeUser(0.01, 100.0);
    hist->GetXaxis()->SetTitle("BCO difference");
    hist->GetYaxis()->SetTitle("Entries");

    if (first_hist) {
      hist->Draw("hist");
      first_hist = false;
    } else {
      hist->Draw("hist same");
    }

    legend->AddEntry(hist, Form("Run %d", run), "l");
    histograms.push_back(hist);
  }

  legend->Draw();

  const TString output_path = Form("%s/h_bco_diff_shifted_%s.pdf",
                                   output_dir.c_str(),
                                   group.label);
  canvas->Print(output_path);
  std::cout << "Saved " << output_path << std::endl;

  for (TH1D* hist : histograms) {
    delete hist;
  }
  delete legend;
  delete canvas;
  return true;
}

}  // namespace

void CompareShiftedHistogramsByGroup(
    const char* input_dir = kDefaultInputDir,
    const char* output_dir = kDefaultOutputDir) {
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(output_dir, true);

  const std::vector<RunGroup> groups = BuildRunGroups();
  const std::vector<int> line_colors = BuildLineColors();

  int num_success = 0;
  for (const RunGroup& group : groups) {
    if (DrawGroup(group, input_dir, output_dir, line_colors)) {
      ++num_success;
    } else {
      std::cerr << "ERROR: Failed to draw " << group.label << std::endl;
    }
  }

  std::cout << "Created " << num_success << " / " << groups.size()
            << " group comparison canvases." << std::endl;
}
