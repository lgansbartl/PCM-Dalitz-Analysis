// Only Plotting

#include "helper.h"

#include <string>
#include <vector>

void plotting()
{

  ///****PRE-SETTINGS *****/
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  // canvas
  std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
  SetCanvasStandardSettings(canvas.get());

  // std::unique_ptr<TCanvas> canvas2pads = std::make_unique<TCanvas>("canvas2pads", "", 800, 800);
  // SetSettingsCanvas_wRatio(canvas2pads.get());

  // header
  std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>();
  SetHeaderSettings(lHeader.get());

  // Open File
  // first ones
  //  TFile* readfile_zo = SafelyOpenRootfile("/Users/lauragans-bartl/master/O2/Plotting/Data/AnalysisResults_LHC23zo.root");
  //  TFile* readfile_zp = SafelyOpenRootfile("/Users/lauragans-bartl/master/O2/Plotting/Data/AnalysisResults_LHC23zp.root");
  //  TFile* readfile_zozp = SafelyOpenRootfile("/Users/lauragans-bartl/master/O2/Plotting/Data/AnalysisResults_LHC23zozp.root");

  // second ones
  TFile* readfile_zo = SafelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults_LHC23zo_0818.root");
  TFile* readfile_zp = SafelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults_LHC23zp_0818.root");
  TFile* readfile_zozp = SafelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults_LHC23zozp_0818.root");

  std::vector<TFile*> vec_readfile = {readfile_zo, readfile_zp, readfile_zozp};
  std::vector<std::string> vec_readfile_name = {"zo", "zp", "zozp"};

  // Getting Files
  std::unique_ptr<THnSparseD> hmpt_same_zo(readfile_zo->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> hmpt_mix_zo(readfile_zo->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision_zo(readfile_zo->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));

  std::unique_ptr<THnSparseD> hmpt_same_zp(readfile_zp->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> hmpt_mix_zp(readfile_zp->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision_zp(readfile_zp->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));

  std::unique_ptr<THnSparseD> hmpt_same_zozp(readfile_zozp->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> hmpt_mix_zozp(readfile_zozp->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision_zozp(readfile_zozp->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));

  std::vector<THnSparseD*> vec_THnsparse_same = {hmpt_same_zo.get(), hmpt_same_zp.get(), hmpt_same_zozp.get()};
  std::vector<THnSparseD*> vec_THnsparse_mix = {hmpt_mix_zo.get(), hmpt_mix_zp.get(), hmpt_mix_zozp.get()};

  // Fittting
  std::unique_ptr<TF1> gaus = std::make_unique<TF1>("gaus", "gaus", 0.47, 0.57);                         // original: 0.40, 0.60
  std::unique_ptr<TF1> fSigAndBack = std::make_unique<TF1>("fSigAndBack", "gaus(0)+ pol2(3)", 0.3, 0.8); // original: 0.40, 0.60 //Zahlen in den Klammern: sagen wo die Paramter fuer die funktionen zu finden sind
  std::unique_ptr<TF1> pol1 = std::make_unique<TF1>("pol1", "pol1", 0., 0.8);
  std::unique_ptr<TF1> pol2 = std::make_unique<TF1>("pol2", "pol2", 0., 0.8);
  std::unique_ptr<TF1> pol3 = std::make_unique<TF1>("pol3", "pol3", 0., 0.8);

  // Vectors for Projection
  std::vector<TH1D*> vec_hProject_same_zo;
  std::vector<TH1D*> vec_hProject_mix_zo;
  std::vector<TH1D*> vec_hProject_same_zp;
  std::vector<TH1D*> vec_hProject_mix_zp;
  std::vector<TH1D*> vec_hProject_same_zozp;
  std::vector<TH1D*> vec_hProject_mix_zozp;

  // Caslculating statistics
  //  double NCollisons_zo = hNCollision_zo->GetBinContent(12);
  //  std::cout << "NCollisons_zo: " << NCollisons_zo << std::endl;
  //  double NCollisons_zp = hNCollision_zp->GetBinContent(12);
  //  std::cout << "NCollisons_zp: " << NCollisons_zp << std::endl;
  //  double NCollisions_cobined_calculated = NCollisons_zo + NCollisons_zp;
  //  std::cout << "NCollisions_cobined_calculated: " << NCollisions_cobined_calculated << std::endl;
  //  double NCollisons_zozp = hNCollision_zozp->GetBinContent(12);
  //  std::cout << "NCollisons_zozp: " << NCollisons_zozp << std::endl;
  //  double NCollisions_total_compared = fabs(NCollisons_zozp - NCollisions_cobined_calculated);
  //  std::cout << "NCollisions_total_compared: " << NCollisions_total_compared << std::endl;

  ////**** THNSPARSE ****////
  for (int NSparse = 0; NSparse < vec_THnsparse_same.size(); NSparse++) {
    // Projections of THnSparse: 0 = minv - histo, 1 = pt - histo
    for (int proj = 0; proj < 2; proj++) { // 0 = minv - histo, 1 = pt - histo
      TH1D* h_same = vec_THnsparse_same[NSparse]->Projection(proj);
      h_same->SetName(Form("h_name_same_%i_%i", NSparse, proj));
      SetHistoStandardSettings1D(h_same);
      SortingVectorN3(NSparse, h_same, vec_hProject_same_zo, vec_hProject_same_zp, vec_hProject_same_zozp);
      h_same->Draw("pe");
      canvas->SaveAs(Form("Plots_%s/hproj_same_%s_%i.png", vec_readfile_name[NSparse].data(), vec_readfile_name[NSparse].data(), proj));

      TH1D* h_mix = vec_THnsparse_mix[NSparse]->Projection(proj);
      h_mix->SetName(Form("h_name_mix_%i_%i", NSparse, proj));
      SetHistoStandardSettings1D(h_mix);
      SortingVectorN3(NSparse, h_same, vec_hProject_mix_zo, vec_hProject_mix_zp, vec_hProject_mix_zozp);
      h_mix->Draw("pe");
      canvas->SaveAs(Form("Plots_%s/hproj_mix_%s_%i.png", vec_readfile_name[NSparse].data(), vec_readfile_name[NSparse].data(), proj));
    }
  }
  std::vector<std::vector<TH1D*>> vec_vecProj_same = {vec_hProject_same_zo, vec_hProject_same_zp, vec_hProject_same_zozp};
  std::vector<std::vector<TH1D*>> vec_vecProj_mix = {vec_hProject_mix_zo, vec_hProject_mix_zp, vec_hProject_mix_zozp};

  // Projections meegamma in pt-slaces
  std::vector<double> vec_pt_min = {0.1, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 2., 3., 5., 18.};
  std::vector<std::string> vec_pt_min_legend = {"0.1", "0.4", "0.5", "0.7", "0.8", "0.9", "1", "2", "3", "5", "18"};
  std::vector<int> vec_rebin = {8, 8, 5, 5, 5, 5, 5, 5, 8, 8};
  double epsilon = 1.e-9; // Hilfsmittel um sicherzustellen, dass der richtige Bin gelesen wird
  std::vector<int> vec_pt_min_bin_same;
  std::vector<int> vec_pt_min_bin_mix;

  // Filling Vectors with bin numbers. Here: vec_hProject_same_zo chosen
  for (int i = 0; i < vec_pt_min.size(); i++) {                        // Filling vector with bins of chosen pt-values
    int j = vec_hProject_same_zo[1]->FindBin(vec_pt_min[i] + epsilon); // epsilon damit die richtigen Bins gelesen werden
    vec_pt_min_bin_same.push_back(j);
    int k = vec_hProject_same_zo[1]->FindBin(vec_pt_min[i] + epsilon);
    vec_pt_min_bin_mix.push_back(k);
  }

  ////****PROJECTIONS ****////
  std::vector<TH1D*> vec_h_same_projm_ptranges_zo;
  std::vector<TH1D*> vec_h_mix_projm_ptranges_zo;
  std::vector<TH1D*> vec_h_same_projm_ptranges_zp;
  std::vector<TH1D*> vec_h_mix_projm_ptranges_zp;
  std::vector<TH1D*> vec_h_same_projm_ptranges_zozp;
  std::vector<TH1D*> vec_h_mix_projm_ptranges_zozp;

  int NRebin = 5;
  // Pt Projections
  for (int NVec = 0; NVec < vec_vecProj_same.size(); NVec++) {
    for (int i = 0; i < vec_pt_min_bin_same.size() - 1; i++) { // Projections of specific pt-ranges
      THnSparseD* hmpt_same_clone = (THnSparseD*)vec_THnsparse_same[NVec]->Clone(Form("hmpt_same_clone_%i", i));
      hmpt_same_clone->GetAxis(1)->SetRange(vec_pt_min_bin_same[i], vec_pt_min_bin_same[i + 1]); // SetRange erwartet Bins. Falls keine Bins SetRangeUsers
      TH1D* h_same = hmpt_same_clone->Projection(0);
      h_same->SetName(Form("h_name_same_projection__%i_%i", NVec, i));
      SetHistoStandardSettings1D(h_same);
      h_same->Rebin(vec_rebin[i]); // kombiniert einzelne Bins zusammen. so viele wie in der Klammer steht.
      SortingVectorN3(NVec, h_same, vec_h_same_projm_ptranges_zo, vec_h_same_projm_ptranges_zp, vec_h_same_projm_ptranges_zozp);
      h_same->Draw("pe");
      canvas->SaveAs(Form("Plots_%s/hproj_same_%f-%f_%i.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1], NRebin));

      // //mixed
      THnSparseD* hmpt_mix_clone = (THnSparseD*)vec_THnsparse_mix[NVec]->Clone(Form("hmpt_mix_clone%i", i));
      hmpt_mix_clone->GetAxis(1)->SetRange(vec_pt_min_bin_mix[i], vec_pt_min_bin_mix[i + 1]);
      TH1D* h_mix = hmpt_mix_clone->Projection(0);
      h_mix->SetName(Form("h_name_mix_projection_%i_%i", NVec, i));
      SetHistoStandardSettings1D(h_mix);
      h_mix->Rebin(vec_rebin[i]);
      SortingVectorN3(NVec, h_mix, vec_h_mix_projm_ptranges_zo, vec_h_mix_projm_ptranges_zp, vec_h_mix_projm_ptranges_zozp);
      h_mix->Draw("pe");
      canvas->SaveAs(Form("Plots_%s/hproj_mix_%f-%f_%i.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1], NRebin));
    }
  }
  std::vector<std::vector<TH1D*>> vec_vecptranges_same = {vec_h_same_projm_ptranges_zo, vec_h_same_projm_ptranges_zp, vec_h_same_projm_ptranges_zozp};
  std::vector<std::vector<TH1D*>> vec_vecptranges_mix = {vec_h_mix_projm_ptranges_zo, vec_h_mix_projm_ptranges_zp, vec_h_mix_projm_ptranges_zozp};

  // Getting Binning of original histos:
  std::vector<TH2D*> vec_h_compared;
  std::vector<double> vecPtBins;
  int NBins;
  NBins = vec_vecptranges_same[0][0]->GetNbinsX();
  vecPtBins.resize(NBins + 1);
  for (int iBin = 1; iBin <= vec_vecptranges_same[0][0]->GetNbinsX() + 1; iBin++) {
    vecPtBins.at(iBin - 1) = vec_vecptranges_same[0][0]->GetBinLowEdge(iBin);
  }

  ////****SIGNAL EXTRACTION****////
  std::vector<TH1D*> vec_h_signal_zo;
  std::vector<TH1D*> vec_h_signal_zp;
  std::vector<TH1D*> vec_h_signal_zozp;
  std::vector<TH1D*> vec_h_mix_scaled_zo;
  std::vector<TH1D*> vec_h_mix_scaled_zp;
  std::vector<TH1D*> vec_h_mix_scaled_zozp;

  for (int NVec = 0; NVec < vec_vecptranges_same.size(); NVec++) {
    for (int i = 0; i < vec_h_same_projm_ptranges_zo.size(); i++) {
      //**1. scaling of mixed
      // 1.1 Defining Histos

      TH1D* h_same = (TH1D*)vec_vecptranges_same[NVec][i]->Clone(Form("h_same_%i_%i", NVec, i));
      TH1D* h_mix = (TH1D*)vec_vecptranges_mix[NVec][i]->Clone(Form("h_mix_%i_%i", NVec, i));

      // 1.2 Scaling Version 1 - one scaling factor
      double skalar_same = h_same->Integral(h_same->FindBin(0.7), h_same->FindBin(0.8));
      double skalar_mix = h_mix->Integral(h_mix->FindBin(0.7), h_mix->FindBin(0.8));
      double scaling_factor = skalar_same / skalar_mix;

      // 1.3 scaling Version 2 - scaling function
      TH1D* h_ratio_same_mix = MakingRatio(h_same, h_mix, "h_ratio_same_mix");
      double binMin = h_ratio_same_mix->FindBin(0.44);
      double binMax = h_ratio_same_mix->FindBin(0.58);
      for (int iBin = binMin; iBin <= binMax; iBin++) {
        h_ratio_same_mix->SetBinContent(iBin, 0);
        h_ratio_same_mix->SetBinError(iBin, 0);
      }
      h_ratio_same_mix->SetAxisRange(0.3, 0.8);
      h_ratio_same_mix->GetYaxis()->SetTitle("same/mix");
      h_ratio_same_mix->Draw("pe");
      h_ratio_same_mix->Fit(pol2.get(), "0M", "", 0.3, 0.8);
      // Plotting Ratio with fit (=function)
      pol2->Draw("SAME");
      SortingVectorN3(NVec, h_ratio_same_mix, vec_h_mix_scaled_zo, vec_h_mix_scaled_zp, vec_h_mix_scaled_zozp);
      canvas->SetLeftMargin(0.09);
      canvas->SaveAs(Form("Plots_%s/h_ratio_same_mix_fit%f-%f.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1]));

      // 1.4 Drawing Scaled mixed Histo
      TH1D* h_mix_scaled = dynamic_cast<TH1D*>(vec_vecptranges_mix[NVec][i]->Clone(Form("h_mix_scaled_%i_%i", NVec, i)));
      SetHistoStandardSettings1D(h_mix_scaled);
      // h_mix_scaled->Scale(scaling_factor);         //sclaing factor
      h_mix_scaled->Multiply(pol2.get()); // scaling function
      h_mix_scaled->Draw("pe");
      SortingVectorN3(NVec, h_mix_scaled, vec_h_mix_scaled_zo, vec_h_mix_scaled_zp, vec_h_mix_scaled_zozp);
      canvas->SaveAs(Form("Plots_%s/h_mix_scaled_%f-%f_new.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1]));

      //**2. Drawing for comparison
      // on one canvas
      std::unique_ptr<TH2D> h_basic_same_comp_scaled_mix = std::make_unique<TH2D>(" h_basic_same_comp_scaled_mix", ";m_{ee#gamma} (GeV/#it{c}^{2}); #it{N}", NBins, 0.3, 0.82, 100, -500, 20e3);
      SetHistoStandardSettings2D(h_basic_same_comp_scaled_mix.get());

      std::unique_ptr<TLegend> legend_same_comp_scaled_mix = std::make_unique<TLegend>(0.26, 0.26, 0.6, 0.34); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
      SetLegendSettings(legend_same_comp_scaled_mix.get());

      std::unique_ptr<TCanvas> canvas2pads = std::make_unique<TCanvas>("canvas2pads", "", 800, 800);
      SetSettingsCanvas_wRatio(canvas2pads.get());
      canvas2pads->cd();

      // Pad 1 (=Compare)
      canvas2pads->cd(1); // wechselt auf das neue Canvas
      h_basic_same_comp_scaled_mix->Draw("");
      h_same->Draw("same");
      h_same->SetLineColor(kBlack);
      h_same->SetMarkerColor(kBlack);
      legend_same_comp_scaled_mix->AddEntry(h_same, "same", "p");
      h_mix_scaled->SetLineColor(kGray + 2);
      h_mix_scaled->SetMarkerColor(kGray + 2);
      h_mix_scaled->SetMarkerStyle(24);
      h_mix_scaled->Draw("same");
      legend_same_comp_scaled_mix->AddEntry(h_mix_scaled, "scaled mixed", "p");

      lHeader->DrawLatexNDC(0.27, 0.50, "ALICE work in progress");
      lHeader->DrawLatexNDC(0.27, 0.45, Form("period: LHC23%s", vec_readfile_name[NVec].data()));
      lHeader->DrawLatexNDC(0.27, 0.40, "pp, #sqrt{#it{s}} = 13.6 TeV, |#it{B}| = 0.2 T");
      lHeader->DrawLatexNDC(0.27, 0.35, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i + 1].data()));
      legend_same_comp_scaled_mix->Draw();
      // canvas2pads->SaveAs(Form("Plots_%s/h_compare_same_scaled_mix_%f-%f.png", vec_readfile_name[NVec].data(),vec_pt_min[i], vec_pt_min[i+1]));

      // Pad 2 (=Ratio)
      canvas2pads->cd(2); // wechselt zurueck
      // Making Ratio of same and scaled mixed
      std::unique_ptr<TH2D> h_basic_signal_ratio = std::make_unique<TH2D>(" h_basic_signal_ratio", ";m_{ee#gamma} (GeV/#it{c}^{2}); same/scaled mixed", NBins, 0.3, 0.82, 100, 0.74, 1.18);
      SetHistoStandardSettings2D(h_basic_signal_ratio.get());

      h_basic_signal_ratio->Draw("");
      // h_basic_signal_ratio->Scale(1., "width");   //hat nichts veraendert
      h_basic_signal_ratio->GetYaxis()->SetNdivisions(505);
      TH1D* h_same_scal_mix_ratio = MakingRatio(h_same, h_mix_scaled, Form("h_K0lPt_ratio_%i_%i", NVec, i));
      h_same_scal_mix_ratio->Draw("same");
      h_same_scal_mix_ratio->GetYaxis()->SetNdivisions(505);
      canvas2pads->SaveAs(Form("Plots_%s/h_ratio_compare_same_scaled_mix%f-%f_new.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1]));
      canvas->cd();
      // 3. Calculating signal
      std::unique_ptr<TLegend> legend_mass = std::make_unique<TLegend>(0.15, 0.65, 0.6, 0.73); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
      SetLegendSettings(legend_mass.get());
      TH1D* h_same_sig = (TH1D*)vec_vecptranges_same[NVec][i]->Clone(Form("h_same_sig_%i_%i", NVec, i));
      h_same_sig->Add(h_same_sig, h_mix_scaled, 1., -1.);
      h_same_sig->GetYaxis()->SetTitle("N");
      SetHistoStandardSettings1D(h_same_sig);
      legend_mass->AddEntry(h_same_sig, Form("LHC23%s", vec_readfile_name[NVec].data()), "p");

      h_same_sig->Fit(gaus.get(), "0M", "", 0.470, 0.570);
      gaus->SetLineColor(kRed);
      // legend_mass->AddEntry(gaus.get(), "Gausian fit component","l");

      fSigAndBack->SetParLimits(1, 0.95 * gaus->GetParameter(1), 1.05 * gaus->GetParameter(1));
      fSigAndBack->SetParLimits(2, 0.5 * gaus->GetParameter(2), 1.05 * gaus->GetParameter(2));
      fSigAndBack->SetLineColor(kAzure);
      fSigAndBack->SetLineWidth(3);
      h_same_sig->Fit(fSigAndBack.get(), "0M", "", 0.2, 0.8);
      legend_mass->AddEntry(fSigAndBack.get(), "pol3 + Gausian fit", "l");
      h_same_sig->GetXaxis()->SetRangeUser(0.28, 0.7);

      // vec_h_signal.push_back(h_same_sig);
      h_same_sig->Draw("pe");

      // gaus->Draw("same");
      lHeader->DrawLatexNDC(0.15, 0.85, "ALICE work in progress");
      lHeader->DrawLatexNDC(0.15, 0.80, "pp, #sqrt{#it{s}} = 13.6 TeV, |#it{B}| = 0.2 T");
      lHeader->DrawLatexNDC(0.15, 0.75, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i + 1].data()));
      fSigAndBack->Draw("same");
      if (i != 0) {
        legend_mass->Draw("");
      }
      SortingVectorN3(NVec, h_same_sig, vec_h_signal_zo, vec_h_signal_zp, vec_h_signal_zozp);
      canvas->SaveAs(Form("Plots_%s/h_signal_%f-%f_new.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1]));

      // 3.1 Comparing Signal with scaled mixed
      std::unique_ptr<TH2D> h_basic_signal_mixed = std::make_unique<TH2D>(" h_basic_signal_mixed", ";m_{ee#gamma} (GeV/#it{c}^{2}); #it{N}", NBins, 0.3, 0.82, 100, -100, 50e3);
      SetHistoStandardSettings2D(h_basic_signal_mixed.get());
      h_basic_signal_mixed->Draw("");
      h_mix_scaled->Draw("same");
      h_same_sig->Draw("same");
      canvas->SaveAs(Form("Plots_%s/h_signal_mixed_%f-%f.png", vec_readfile_name[NVec].data(), vec_pt_min[i], vec_pt_min[i + 1]));
    }
  }

  // 3.1 Only at Low pt
  std::unique_ptr<TLegend> legend_mass_lowfield = std::make_unique<TLegend>(0.14, 0.15, 0.6, 0.23); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
  SetLegendSettings(legend_mass_lowfield.get());
  legend_mass_lowfield->AddEntry(vec_h_signal_zo[0], "LHC23zo", "p");
  fSigAndBack->SetParLimits(1, 0.95 * gaus->GetParameter(1), 1.05 * gaus->GetParameter(1));
  fSigAndBack->SetParLimits(2, 0.5 * gaus->GetParameter(2), 1.05 * gaus->GetParameter(2));
  fSigAndBack->SetLineColor(kAzure);
  fSigAndBack->SetLineWidth(3);
  vec_h_signal_zo[0]->Fit(fSigAndBack.get(), "0M", "", 0.2, 0.8);
  legend_mass_lowfield->AddEntry(fSigAndBack.get(), "pol3 + Gausian fit", "l");
  vec_h_signal_zo[0]->Draw("pe");
  fSigAndBack->Draw("same");
  lHeader->DrawLatexNDC(0.55, 0.35, "ALICE work in progress");
  lHeader->DrawLatexNDC(0.15, 0.30, "pp, #sqrt{#it{c}} = 13.6 TeV, |#it{B}| = 0.2 T");
  lHeader->DrawLatexNDC(0.15, 0.25, "0.1 GeV/#it{c} < #it{p}_{T} < 0.4 GeV/#it{c}");
  legend_mass_lowfield->Draw("");
  canvas->SaveAs("Plots_zo/h_signal_lowfield.png");

  // 4 Drawing all Datasets in one canvas
  // low field
  double maxYAxis = 3200;
  std::unique_ptr<TH2D> h_basic_signal_compared_0 = std::make_unique<TH2D>(" h_basic_signal_compared", ";m_{ee#gamma} (GeV/#it{c}^{2}); #it{N}", NBins, 0.3, 0.82, 100, -500, maxYAxis);
  SetHistoStandardSettings2D(h_basic_signal_compared_0.get());

  std::unique_ptr<TLegend> legend_signal_0 = std::make_unique<TLegend>(0.26, 0.66, 0.6, 0.78); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
  SetLegendSettings(legend_signal_0.get());

  h_basic_signal_compared_0->Draw("");
  vec_h_signal_zo[0]->SetLineColor(kAzure);
  vec_h_signal_zo[0]->SetMarkerColor(kAzure);
  vec_h_signal_zo[0]->Draw("same");
  legend_signal_0->AddEntry(vec_h_signal_zo[0], "LHC23zo", "p");
  vec_h_signal_zp[0]->SetLineColor(kGreen + 2);
  vec_h_signal_zp[0]->SetMarkerColor(kGreen + 2);
  vec_h_signal_zp[0]->Draw("same");
  legend_signal_0->AddEntry(vec_h_signal_zp[0], "LHC23zp", "p");
  vec_h_signal_zozp[0]->SetLineColor(1);
  vec_h_signal_zozp[0]->SetMarkerColor(1);
  vec_h_signal_zozp[0]->Draw("same");
  legend_signal_0->AddEntry(vec_h_signal_zozp[0], "combined", "p");

  lHeader->DrawLatexNDC(0.27, 0.90, "ALICE work in progress");
  lHeader->DrawLatexNDC(0.27, 0.85, "pp, #sqrt{#it{s}} = 13.6 TeV, |#it{B}| = 0.2 T");
  lHeader->DrawLatexNDC(0.27, 0.80, "0.1 GeV/#it{c} < #it{p}_{T} < 0.4 GeV/#it{c}");

  legend_signal_0->Draw("");
  canvas->SaveAs("Plots_general/h_signal_compared_lowfield.png");

  // rest
  for (int NElements = 1; NElements < vec_h_signal_zo.size(); NElements++) {
    double maxYAxis = 3000;
    std::unique_ptr<TH2D> h_basic_signal_compared = std::make_unique<TH2D>(" h_basic_signal_compared", ";m_{ee#gamma} (GeV/#it{c}^{2}); #it{N}", NBins, 0.3, 0.82, 100, -500, maxYAxis);
    SetHistoStandardSettings2D(h_basic_signal_compared.get());

    std::unique_ptr<TLegend> legend_signal = std::make_unique<TLegend>(0.39, 0.63, 0.6, 0.75); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    SetLegendSettings(legend_signal.get());

    h_basic_signal_compared->Draw("");
    vec_h_signal_zo[NElements]->SetLineColor(kAzure);
    vec_h_signal_zo[NElements]->SetMarkerColor(kAzure);
    vec_h_signal_zo[NElements]->Draw("same");
    legend_signal->AddEntry(vec_h_signal_zo[NElements], "LHC23zo", "p");
    vec_h_signal_zp[NElements]->SetLineColor(kGreen + 2);
    vec_h_signal_zp[NElements]->SetMarkerColor(kGreen + 2);
    vec_h_signal_zp[NElements]->Draw("same");
    legend_signal->AddEntry(vec_h_signal_zp[NElements], "LHC23zp", "p");
    vec_h_signal_zozp[NElements]->SetLineColor(1);
    vec_h_signal_zozp[NElements]->SetMarkerColor(1);
    vec_h_signal_zozp[NElements]->Draw("same");
    legend_signal->AddEntry(vec_h_signal_zozp[NElements], "combined", "p");

    lHeader->DrawLatexNDC(0.40, 0.88, "ALICE work in progress");
    lHeader->DrawLatexNDC(0.40, 0.83, "pp, #sqrt{#it{s}} = 13.6 TeV, |#it{B}| = 0.2 T");
    lHeader->DrawLatexNDC(0.40, 0.78, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[NElements].data(), vec_pt_min_legend[NElements + 1].data()));

    legend_signal->Draw("");
    canvas->SaveAs(Form("Plots_general/h_signal_compared_%f-%f_sfactor.png", vec_pt_min[NElements], vec_pt_min[NElements + 1]));
    // canvas->SaveAs(Form("Plots_general/h_signal_compared_%f-%f_sfunction.png", vec_pt_min[NElements], vec_pt_min[NElements+1]));
    vec_h_compared.push_back(h_basic_signal_compared.get()); // irgendetwas geht schief!!!!
  }
}

int main(void)
{

  plotting();

  return 0;
}
