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

#include "Utility/sPhenixStyle.C"

void PutCaption(double first_margin = -0.1, double pos_x = 0.2, double pos_y = 0.9);

namespace {

struct MultiplicityConfig {
  const char* dir_suffix;
  const char* legend_label;
  int line_color;
};

const char* kDefaultBaseOutputDir =
    "/sphenix/tg/tg01/commissioning/INTT/work/ryotaro/TimingResolution/output";

std::vector<MultiplicityConfig> BuildMultiplicityConfigs() {
  return {
      // {"mul10to100", "10 < trigger timing hits < 100", kRed + 1},
      // {"mul100to200", "100 < trigger timing hits < 200", kBlue + 1},
      // {"mul200to1000", "200 < trigger timing hits < 1000", kGreen + 2},
      {"mul10to100", "Low multiplicity", kRed + 1},
      {"mul100to200", "Middle multiplicity", kBlue + 1},
      {"mul200to1000", "High multiplicity", kGreen + 2},
  };
}

TH1D* LoadSolvedHistogram(const std::string& base_output_dir,
                          const std::string& group_label,
                          const MultiplicityConfig& multiplicity) {
  const TString input_path = Form(
      "%s/shifted_sum_reconstruction_%s_%s/solved_n.root",
      base_output_dir.c_str(),
      group_label.c_str(),
      multiplicity.dir_suffix);

  TFile* input_file = TFile::Open(input_path, "READ");
  if (!input_file || input_file->IsZombie()) {
    std::cerr << "ERROR: Could not open " << input_path << std::endl;
    delete input_file;
    return nullptr;
  }

  TH1D* input_hist = dynamic_cast<TH1D*>(input_file->Get("h_solved_n"));
  if (input_hist == nullptr) {
    std::cerr << "ERROR: h_solved_n is missing in " << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(
      input_hist->Clone(Form("h_solved_n_%s_%s", group_label.c_str(), multiplicity.dir_suffix)));
  if (hist == nullptr) {
    std::cerr << "ERROR: Failed to clone h_solved_n from " << input_path << std::endl;
    input_file->Close();
    delete input_file;
    return nullptr;
  }

  hist->SetDirectory(nullptr);
  input_file->Close();
  delete input_file;
  return hist;
}

bool DrawGroup(const std::string& group_label,
               const std::string& base_output_dir,
               const std::vector<MultiplicityConfig>& multiplicities) {
  TCanvas* canvas =
      new TCanvas(Form("c_compare_solved_n_%s", group_label.c_str()), group_label.c_str(), 900, 700);
  // canvas->SetMargin(0.12, 0.04, 0.11, 0.06);
  // canvas->SetTicks(1, 1);

  // TLegend* legend = new TLegend(0.52, 0.68, 0.88, 0.88);
  // TLegend* legend = new TLegend(0.4, 0.6, 0.65, 0.8);
  // TLegend* legend = new TLegend(0.2, 0.5, 0.45, 0.7);
  TLegend* legend = new TLegend(0.6, 0.7, 0.8, 0.94);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);
  // legend->SetTextSize(0.03);

  std::vector<TH1D*> histograms;
  histograms.reserve(multiplicities.size());

  bool first_hist = true;
  for (std::size_t index = 0; index < multiplicities.size(); ++index) {
    const MultiplicityConfig& multiplicity = multiplicities[index];
    TH1D* hist = LoadSolvedHistogram(base_output_dir, group_label, multiplicity);
    if (hist == nullptr) {
      for (TH1D* loaded_hist : histograms) {
        delete loaded_hist;
      }
      delete legend;
      delete canvas;
      return false;
    }

    // hist->SetStats(0);
    // hist->SetTitle("");
    hist->SetLineColor(multiplicity.line_color);
    // hist->SetLineWidth(2);
    hist->GetXaxis()->SetTitle("hit timing [1/6 BCO]");
    hist->GetYaxis()->SetTitle("Number of INTT hits (event averaged)");
    hist->GetYaxis()->SetRangeUser(-5.0, 300.0);

    if (first_hist) {
      hist->Draw("hist");
      first_hist = false;
    } else {
      hist->Draw("hist same");
    }

    legend->AddEntry(hist, multiplicity.legend_label, "l");
    PutCaption(-0.1,0.2,0.94);
    histograms.push_back(hist);
  }

  legend->Draw();

  const TString output_path = Form(
      "%s/fine_hit_timing_distribution_multiplicity_comparison_%s.pdf",
      base_output_dir.c_str(),
      group_label.c_str());
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

void CompareSolvedNByMultiplicity(const char* base_output_dir = kDefaultBaseOutputDir) {
  SetsPhenixStyle();
  
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(base_output_dir, true);

  const std::vector<MultiplicityConfig> multiplicities = BuildMultiplicityConfigs();
  const std::vector<std::string> groups = {"Group1", "Group1_wo43280", "Group2"};

  int num_success = 0;
  for (const std::string& group_label : groups) {
    if (DrawGroup(group_label, base_output_dir, multiplicities)) {
      ++num_success;
    } else {
      std::cerr << "ERROR: Failed to draw " << group_label << std::endl;
    }
  }

  std::cout << "Created " << num_success << " / " << groups.size()
            << " multiplicity comparison canvases." << std::endl;
}


void PutCaption(double first_margin = -0.1, double pos_x = 0.2, double pos_y = 0.9){
/*from Genki*/
  TLatex* tex = new TLatex();
  double line_height = 0.055; // no need to change
  // double first_margin = -0.1; // the margin from the top of the canvas to the top line (date). You can modify
  // double pos_x = 0.2; // you can modify
  // double pos_y = 1.0  + first_margin; // - line_height;

  /* Date, you may want to change 0.65 to pos_x */
  //tex->DrawLatexNDC( pos_x, pos_y,
  //           string("#it{" + GetDate() + "}").c_str() );

  // sPHENIX Internal or sPHENIX Prelimnary
  pos_y -= line_height;
//   tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Preliminary" );
  // tex->DrawLatexNDC( pos_x, 0.955, "#it{#bf{sPHENIX}} Internal" );
  // tex->DrawLatexNDC( pos_x, 0.955, "#it{#bf{sPHENIX}} Preliminary" );
  //tex->DrawLatexNDC( pos_x, pos_y, string( string("#it{#bf{sPHENIX}} Preliminary ") + "#it{" + GetDate() + "}").c_str() );
  //tex->DrawLatexNDC( pos_x, pos_y, string("#it{#bf{sPHENIX}} Preliminary #it{05/16/2025}").c_str() );

  // tex->DrawLatexNDC( 0.77, 0.955, string("#it{05/13/2026}").c_str() );
  // tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Internal #it{05/13/2026}" );
  tex->DrawLatexNDC( pos_x, pos_y, "#it{05/13/2026}" );
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, "#it{#bf{sPHENIX}} Internal" );
  
  // p+p 200 GeV
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, ("Run-24 #it{p+p}") );
  pos_y -= line_height ;// + 0.03;
  tex->DrawLatexNDC( pos_x, pos_y, ("#sqrt{s} = 200 GeV") );
  // tex->DrawLatexNDC( pos_x, 0.88, ("#it{p+p} #sqrt{s} = 200 GeV") );


}

