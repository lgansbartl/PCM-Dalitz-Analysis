// Plotting fuer FSP-Meeting: etaSignal zeigen von Pion und Eta

#include "helper.h"

void plotting_nf_data()
{

  ///****PRE-SETTINGS *****/
  gStyle->SetOptStat(1);
  TGaxis::SetMaxDigits(3);

  // canvas
  std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
  SetCanvasStandardSettings(canvas.get());

  // header
  std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>();
  SetHeaderSettings(lHeader.get());

  // TFile* usedFile;

  // open file
  //  TFile* readfile = SafelyOpenRootfile("/Users/lauragans-bartl/master/data/AnalysisResults_nf_dd.root");
  TFile* readfile = SafelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC25ddsmall.root");

  std::unique_ptr<THnSparseD> hmpt_same_zozp(readfile->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> hmpt_mix_zozp(readfile->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision_zozp(readfile->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));

  // getting file

  // Defining vectors for projection
  std::vector<TH1D*> vec_hProject_same_zozp;
  std::vector<TH1D*> vec_hProject_mix_zozp;
  std::vector<std::string> vec_proj_type = {"minv", "pt"};

  std::vector<double> vec_pt_min = {0.1, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 2., 3., 5., 18.};

  // Fittting
  std::unique_ptr<TF1> gaus = std::make_unique<TF1>("gaus", "gaus", 0.47, 0.57);                           // original: 0.40, 0.60
  std::unique_ptr<TF1> gaus2 = std::make_unique<TF1>("gaus2", "gaus", 0.05, 0.17);                         // original: 0.40, 0.60
  std::unique_ptr<TF1> fSigAndBack = std::make_unique<TF1>("fSigAndBack", "gaus(0)+ pol2(3)", 0.40, 0.63); // original: 0.40, 0.60 //Zahlen in den Klammern: sagen wo die Paramter fuer die funktionen zu finden sind
  std::unique_ptr<TF1> pol2 = std::make_unique<TF1>("pol2", "pol2", 0., 0.8);

  // Massepeak in Pt
  std::unique_ptr<TH1D> h_MmaxPTEta = std::make_unique<TH1D>(" h_MmaxPTEta", ";#it{p}_{T} (GeV/#it{c}); #LT #it{m}_{#eta} #GT (GeV/#it{c}^{2})", vec_pt_min.size() - 1, vec_pt_min.data());
  SetHistoStandardSettings1D(h_MmaxPTEta.get(), 1.2, 1.6);

  std::unique_ptr<TH1D> h_MmaxPTPion = std::make_unique<TH1D>(" h_MmaxPTPion", ";#it{p}_{T} (GeV/#it{c}); #LT #it{m}_{#pi_{0}} #GT (GeV/#it{c}^{2})", vec_pt_min.size() - 1, vec_pt_min.data());
  SetHistoStandardSettings1D(h_MmaxPTPion.get(), 1.2, 1.7);

  // Projections aus THNSparse
  for (int proj = 0; proj < 2; proj++) { // 0 = minv - histo, 1 = pt - histo

    TH1D* h_same = hmpt_same_zozp->Projection(proj);
    h_same->SetName(Form("h_name_same_%i", proj));
    SetHistoStandardSettings1D(h_same);
    h_same->Draw("pe");
    vec_hProject_same_zozp.push_back(h_same);
    // canvas->SaveAs(Form("Plots_nf_data/hproj_same_%s.png", vec_proj_type[proj].data()));
    TH1D* h_mix = hmpt_mix_zozp->Projection(proj);
    h_mix->SetName(Form("h_name_mix_%i", proj));
    SetHistoStandardSettings1D(h_mix);
    h_mix->Draw("pe");
    vec_hProject_mix_zozp.push_back(h_mix);

    // canvas->SaveAs(Form("Plots_nf_data/hproj_mix_%s.png", vec_proj_type[proj].data()));
  }

  //****PT-Ptojections****//
  // Pt-Bins werden gesetzt
  std::vector<std::string> vec_pt_min_legend = {"0.1", "0.4", "0.5", "0.7", "0.8", "0.9", "1", "2", "3", "5", "18"};
  std::vector<int> vec_rebin = {8, 8, 5, 5, 5, 5, 5, 5, 8, 8};
  double epsilon = 1.e-9; // Hilfsmittel um sicherzustellen, dass der richtige Bin gelesen wird
  std::vector<int> vec_pt_min_bin_same;
  std::vector<int> vec_pt_min_bin_mix;
  // Passende Bin Nummern werden aus dem pt HIstogramm geholt
  for (int i = 0; i < vec_pt_min.size(); i++) {                          // Filling vector with bins of chosen pt-values
    int j = vec_hProject_same_zozp[1]->FindBin(vec_pt_min[i] + epsilon); // epsilon damit die richtigen Bins gelesen werden
    vec_pt_min_bin_same.push_back(j);
    int k = vec_hProject_mix_zozp[1]->FindBin(vec_pt_min[i] + epsilon);
    vec_pt_min_bin_mix.push_back(k);
  }

  int NRebin = 5;
  // projection der MAsseverteilung mit den einzelnen pt Bereichen
  std::vector<TH1D*> vec_h_same_projm_ptranges_zozp;
  std::vector<TH1D*> vec_h_mix_projm_ptranges_zozp;
  for (int i = 0; i < vec_pt_min_bin_same.size() - 1; i++) { // Projections of specific pt-ranges
    // same
    THnSparseD* hmpt_same_clone = (THnSparseD*)hmpt_same_zozp->Clone(Form("hmpt_same_clone_%i", i));
    hmpt_same_clone->GetAxis(1)->SetRange(vec_pt_min_bin_same[i], vec_pt_min_bin_same[i + 1]); // SetRange erwartet Bins. Falls keine Bins SetRangeUsers
    TH1D* h_same = hmpt_same_clone->Projection(0);                                             // masseverteilung wird rausgeholt
    h_same->SetName(Form("h_name_same_projection_%i", i));
    SetHistoStandardSettings1D(h_same);
    h_same->Rebin(vec_rebin[i]); // kombiniert einzelne Bins zusammen. so viele wie in der Klammer steht.
    h_same->Draw("pe");
    vec_h_same_projm_ptranges_zozp.push_back(h_same);
    // canvas->SaveAs(Form("Plots_nf_data/hproj_same_%0.1f-%0.1f_Rebin%i.png", vec_pt_min[i], vec_pt_min[i+1], NRebin));

    // mixed
    THnSparseD* hmpt_mix_clone = (THnSparseD*)hmpt_mix_zozp->Clone(Form("hmpt_mix_clone%i", i));
    hmpt_mix_clone->GetAxis(1)->SetRange(vec_pt_min_bin_mix[i], vec_pt_min_bin_mix[i + 1]);
    TH1D* h_mix = hmpt_mix_clone->Projection(0);
    h_mix->SetName(Form("h_name_mix_projection_%i", i));
    SetHistoStandardSettings1D(h_mix);
    h_mix->Rebin(vec_rebin[i]);
    h_mix->Draw("pe");
    vec_h_mix_projm_ptranges_zozp.push_back(h_mix);
    // canvas->SaveAs(Form("Plots_nf_data/hproj_mix_%0.1f-%0.1f_Rebin%i.png", vec_pt_min[i], vec_pt_min[i+1], NRebin));
  }

  // Getting Binning of original histos:
  std::vector<TH2D*> vec_h_compared;
  std::vector<double> vecPtBins;
  int NBins;
  NBins = vec_h_same_projm_ptranges_zozp[0]->GetNbinsX();
  vecPtBins.resize(NBins + 1);
  for (int iBin = 1; iBin <= vec_h_same_projm_ptranges_zozp[0]->GetNbinsX() + 1; iBin++) {
    vecPtBins.at(iBin - 1) = vec_h_same_projm_ptranges_zozp[0]->GetBinLowEdge(iBin);
  }

  // Defining Vectors for next loop
  std::vector<double> VecEtaSoverBMixed;
  std::vector<double> VecPionSoverBMixed;
  //****Extraction****//
  for (int i = 0; i < vec_h_same_projm_ptranges_zozp.size(); i++) {
    //**1.1 Defining Histos
    TH1D* h_same = (TH1D*)vec_h_same_projm_ptranges_zozp[i]->Clone(Form("h_same_%i", i));
    TH1D* h_mix = (TH1D*)vec_h_mix_projm_ptranges_zozp[i]->Clone(Form("h_mix_%i", i));

    //**1.2 Scaling
    //**1.2.1 Scaling Version 1 - one scaling factor
    double skalar_same = h_same->Integral(h_same->FindBin(0.7), h_same->FindBin(0.8));
    double skalar_mix = h_mix->Integral(h_mix->FindBin(0.7), h_mix->FindBin(0.8));
    double scaling_factor = skalar_same / skalar_mix;

    //**1.2.2 scaling Version 2 - scaling function
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
    //**1.3 Drawing Scaled mixed Histo
    TH1D* h_mix_scaled = (TH1D*)vec_h_mix_projm_ptranges_zozp[i]->Clone(Form("h_mix_scaled_%i", i));
    SetHistoStandardSettings1D(h_mix_scaled);
    // h_mix_scaled->Scale(scaling_factor);         //scaling factor
    h_mix_scaled->Multiply(pol2.get()); // scaling function
    h_mix_scaled->Draw("pe");
    // canvas->SaveAs(Form("Plots_nf_data/h_mix_scaled_%0.1f-%0.1f.png",vec_pt_min[i], vec_pt_min[i+1]));

    TH1* hSameScaledMixRatio = MakingRatio(h_same, h_mix_scaled, "hSameScaledMixRatio");
    hSameScaledMixRatio->Draw("");
    canvas->SaveAs(Form("Plots_nf_data/hSameScaledMixRatio%0.1f-%0.1f.png", vec_pt_min[i], vec_pt_min[i + 1]));

    // TH1* hSameMixRatio = MakingRatio(h_same, h_mix, "hSameMixRatio");
    // hSameMixRatio->Draw("");
    // canvas->SaveAs(Form("Plots_nf_data/hSameMixRatio_%0.1f-%0.1f.png",vec_pt_min[i], vec_pt_min[i+1]));

    //** 2. Calculating Eta-Signal
    std::unique_ptr<TLegend> legend_mass = std::make_unique<TLegend>(0.4, 0.35, 0.70, 0.51); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    SetLegendSettings(legend_mass.get());
    TH1D* h_signal = (TH1D*)vec_h_same_projm_ptranges_zozp[i]->Clone(Form("h_same_sig_%i", i));
    h_signal->Add(h_signal, h_mix_scaled, 1., -1.);
    h_signal->GetYaxis()->SetTitle("#it{N}");
    h_signal->GetXaxis()->SetTitle("#it{m}_{ee#gamma} (GeV/#it{c}^{2})");
    SetHistoStandardSettings1D(h_signal);

    h_signal->Fit(gaus.get(), "0M", "", 0.470, 0.570);
    gaus->SetLineColor(kRed);
    //**Fitting Eta
    fSigAndBack->SetParLimits(1, 0.95 * gaus->GetParameter(1), 1.05 * gaus->GetParameter(1));
    fSigAndBack->SetParLimits(2, 0.5 * gaus->GetParameter(2), 1.05 * gaus->GetParameter(2));
    fSigAndBack->SetLineColor(kRed + 1);
    fSigAndBack->SetLineWidth(3);
    fSigAndBack->SetNpx(10000);
    h_signal->Fit(fSigAndBack.get(), "0M", "", 0.2, 0.8);
    // h_signal->GetXaxis()->SetRangeUser(0.28, 0.7);

    //**Einschub: Massepeaks werden rausgezogen, em die Abweichung von der erwareten Masse zu plotten. Plot: nach dem Loop
    double m_max_eta = gaus->GetParameter(1);
    double m_max_eta_error = gaus->GetParError(1);
    h_MmaxPTEta->SetBinContent(i + 1, m_max_eta);
    h_MmaxPTEta->SetBinError(i + 1, m_max_eta_error);

    //**Fitting Pion with Steffis function
    TF1* FitPion = FitPi0Shape_Save(h_signal, 0.02, 0.25);
    FitPion->SetLineColor(kBlue);
    FitPion->SetLineWidth(3);
    FitPion->Draw("same");

    //**Einschub: Massepeaks werden rausgezogen, em die Abweichung von der erwareten Masse zu plotten. Plot: nach dem Loop
    double m_max_pion = FitPion->GetParameter(4);
    double m_max_pion_error = FitPion->GetParError(4);

    h_MmaxPTPion->SetBinContent(i + 1, m_max_pion);
    h_MmaxPTPion->SetBinError(i + 1, m_max_pion_error);

    h_signal->Draw("pe");
    fSigAndBack->Draw("same");

    // Linien zum Einzeichnen vom Massenbereich
    //  std::unique_ptr<TLine> linePionMin = std::make_unique<TLine> (0.06, -500, 0.06, 2000);
    //  linePionMin->SetLineWidth(3);
    //  linePionMin->SetLineColor(kAzure-4);
    //  linePionMin->Draw("same");

    // std::unique_ptr<TLine> linePionMax = std::make_unique<TLine> (0.18, -500, 0.18, 2000);
    // linePionMax->SetLineWidth(3);
    // linePionMax->SetLineColor(kAzure-4);
    // linePionMax->Draw("same");

    // std::unique_ptr<TLine> lineEtaMin = std::make_unique<TLine> (0.46, -500, 0.46, 2000);
    // lineEtaMin->SetLineWidth(3);
    // lineEtaMin->SetLineColor(kRed-8);
    // lineEtaMin->Draw("same");

    // std::unique_ptr<TLine> lineEtaMax = std::make_unique<TLine> (0.58, -500, 0.58, 2000);
    // lineEtaMax->SetLineWidth(3);
    // lineEtaMax->SetLineColor(kRed-8);
    // lineEtaMax->Draw("same");

    // h_signal->GetYaxis()->SetRangeUser(-100., 50e3);
    h_same->Draw("same");
    h_same->SetMarkerColor(kGray);
    h_same->SetLineColor(kGray);
    legend_mass->AddEntry(h_signal, "same event - mixed event", "p");
    legend_mass->AddEntry(h_same, "same event", "p");
    legend_mass->AddEntry(FitPion, "Crystall ball fit", "l");
    legend_mass->AddEntry(fSigAndBack.get(), "pol2 + Gausian fit", "l");
    // DrawingHeaderStandardLines(0.4, 0.70, 0.65, 0.60);
    // lHeader->DrawLatexNDC(0.4, 0.55, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i+1].data()));

    // legend_mass->Draw("same");
    canvas->SaveAs(Form("Plots_nf_data/h_signal_eta_pion%0.1f-%0.1f.png", vec_pt_min[i], vec_pt_min[i + 1]));
    canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/h_signal_eta_pion%0.1f-%0.1f.pdf", vec_pt_min[i], vec_pt_min[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/h_signal_eta_pion_massspectrum%0.1f-%0.1f.pdf", vec_pt_min[i], vec_pt_min[i+1]));        //fuer den Mas

    //**NUR Pion Signal
    h_signal->GetXaxis()->SetRangeUser(0., 0.2);
    // Legende fuer 0.8< pt < 0.9
    //  lHeader->DrawLatexNDC(0.16,0.70, "ALICE work in progress");
    //  lHeader->DrawLatexNDC(0.16,0.65, "LHC23zo,zp");
    //  lHeader->DrawLatexNDC(0.16, 0.60, "pp, #sqrt{#it{s}} = 13.6 TeV, B = 0.2 T");
    //  lHeader->DrawLatexNDC(0.16,0.55, "#pi^{0} #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");
    //  lHeader->DrawLatexNDC(0.16, 0.50, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i+1].data()));

    canvas->SaveAs(Form("Plots_nf_data/hSignalPion_%0.1f_%0.1f.png", vec_pt_min[i], vec_pt_min[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/hSignalPion_%0.1f-%0.1f.pdf", vec_pt_min[i], vec_pt_min[i+1]));

    h_signal->GetXaxis()->SetRangeUser(0.3, 0.7);
    h_signal->GetYaxis()->SetRangeUser(-500, 550); // fuer den kleinen Eta Bereich

    //**NUR ETA Signal
    // Legende fuer 1< pt < 2
    // lHeader->DrawLatexNDC(0.13,0.90, "ALICE work in progress");
    // lHeader->DrawLatexNDC(0.13,0.85, "LHC23zo,zp");
    // lHeader->DrawLatexNDC(0.13, 0.80, "pp, #sqrt{#it{s}} = 13.6 TeV, B = 0.2 T");
    // lHeader->DrawLatexNDC(0.13,0.75, "#eta #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");
    // lHeader->DrawLatexNDC(0.13, 0.70, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i+1].data()));

    // Legende fuer 0.8< pt < 0.9
    DrawingHeaderStandardLines(0.13, 0.91, 0.28, 0.23);
    lHeader->DrawLatexNDC(0.13, 0.18, "#eta #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");
    lHeader->DrawLatexNDC(0.13, 0.13, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vec_pt_min_legend[i].data(), vec_pt_min_legend[i + 1].data()));

    canvas->SaveAs(Form("Plots_nf_data/hSignalEta_%0.1f_%0.1f.png", vec_pt_min[i], vec_pt_min[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/hSignalEta_%0.1f-%0.1f.pdf", vec_pt_min[i], vec_pt_min[i+1]));

    //**3.Signal/Background (Wert fuer den jeweiligen pt-Bereich)
    //**3.1: fuer Eta
    double etaSignal = h_signal->Integral(h_signal->FindBin(0.5), h_signal->FindBin(0.55)); //->FindBin(0.46), h_signal->FindBin(0.57));
    //**Version1: scaled mixed als Hintergrund
    double EtaBackgroundMixed = h_mix_scaled->Integral(h_mix_scaled->FindBin(0.5), h_mix_scaled->FindBin(0.55));
    double SoverB_mixed_eta = etaSignal / EtaBackgroundMixed;
    VecEtaSoverBMixed.push_back(SoverB_mixed_eta);
    //**Version2: same Signal als Hintergrund
    // h_same->Add(h_same, h_signal, 1, -1);
    // double EtaBackgroundSameSignal = h_same->Integral(h_same->FindBin(0.5), h_same->FindBin(0.55));
    // double EtaSoverB_same_signal = etaSignal/EtaBackgroundSameSignal;

    //**3.2 fuer Pion
    double pionSignal = h_signal->Integral(h_signal->FindBin(0.08), h_signal->FindBin(0.14));
    //**Version1: scaled mixed als Huntergrund
    double PionBackgroundMixed = h_mix_scaled->Integral(h_mix_scaled->FindBin(0.08), h_mix_scaled->FindBin(0.14));
    double pionSoverB_mixed = pionSignal / PionBackgroundMixed;
    VecPionSoverBMixed.push_back(pionSoverB_mixed);
    //**Version2: same-etaSignal als Hintergrund
    // h_same->Add(h_same, h_signal, 1, -1);
    // double pionBackgroundSameSignal = h_same->Integral(h_same->FindBin(0.6), h_same->FindBin(0.18));
    // double PionSoverB_same_signal = pionSignal/pionBackgroundSameSignal;
  }

  //**4.4 Signal/Background
  //**4.4.1 Eta
  canvas->SetLogy(1);
  std::unique_ptr<TH1D> h_EtaSoverB = std::make_unique<TH1D>("h_EtaSoverB", " ;#it{p}_{T} (GeV/#it{c}); S/B", vec_pt_min.size() - 1, vec_pt_min.data());
  SetHistoStandardSettings1D(h_EtaSoverB.get());
  h_EtaSoverB->GetXaxis()->SetRangeUser(0.8, 6);
  h_EtaSoverB->GetYaxis()->SetRangeUser(1e-2, 12);
  for (int i = 0; i < vec_pt_min.size(); i++) {
    if (vec_pt_min[i] > 0.8) {
      h_EtaSoverB->SetBinContent(h_EtaSoverB->GetXaxis()->FindBin(vec_pt_min[i]), VecEtaSoverBMixed[i]);
      h_EtaSoverB->SetBinError(h_EtaSoverB->GetXaxis()->FindBin(vec_pt_min[i]), 10e-10);
    }
  }
  h_EtaSoverB->SetMarkerColor(kRed + 1);
  h_EtaSoverB->SetLineColor(kRed + 1);

  h_EtaSoverB->Draw("");

  // canvas->SaveAs("Plots_nf_data/h_EtaSoverB.png");
  // canvas->SaveAs("/Users/lauragans-bartl/Downloads/h_EtaSoverB.pdf");

  //**4.4.2 Pion
  std::unique_ptr<TH1D> h_PionSoverB = std::make_unique<TH1D>("h_PionSoverB", " ;#it{p}_{T} (GeV/#it{c}); S/B", vec_pt_min.size() - 1, vec_pt_min.data());
  SetHistoStandardSettings1D(h_PionSoverB.get());
  for (int i = 0; i < vec_pt_min.size() - 1; i++) {
    h_PionSoverB->SetBinContent(i + 1, VecPionSoverBMixed[i]);
    h_PionSoverB->SetBinError(i + 1, 10e-10);
  }
  h_PionSoverB->SetMarkerColor(kBlue);
  h_PionSoverB->SetLineColor(kBlue);
  h_PionSoverB->Draw("same");

  // canvas->SaveAs("Plots_nf_data/h_PionSoverB.png");

  //**Beide auf einem Canvas
  std::unique_ptr<TLegend> legendSoverB = std::make_unique<TLegend>(0.36, 0.24, 0.65, 0.32); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
  SetLegendSettings(legendSoverB.get());
  DrawingHeaderStandardLines(0.4, 0.50, 0.45, 0.40);
  lHeader->DrawLatexNDC(0.4, 0.35, "#pi^{0}, #eta, #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");
  legendSoverB->AddEntry(h_PionSoverB.get(), "#pi^{0}: 0.06 GeV/#it{c}^{2} < #it{m}_{ee#gamma} < 0.18 GeV/#it{c}^{2}", "p");
  legendSoverB->AddEntry(h_EtaSoverB.get(), "#eta : 0.5 GeV/#it{c}^{2} < #it{m}_{ee#gamma} < 0.55 GeV/#it{c}^{2}", "p");

  // h_PionSoverB_signal->Draw("same");
  legendSoverB->Draw("same");
  // canvas->SaveAs("Plots_nf_data/h_PionSoverB.png");
  // SavingPdfPng("h_PionEtaSoverB");
  canvas->SaveAs("Plots_nf_data/h_PionEtaSoverB.png");
  // canvas->SaveAs("/Users/lauragans-bartl/Downloads/h_PionEtaSoverB.pdf");
  // canvas->SaveAs("/Users/lauragans-bartl/Downloads/h_PionSoverB.pdf");

  canvas->SetLogy(0);

  //***Massverschiebung
  //**Eta
  canvas->SetLeftMargin(0.11);
  std::unique_ptr<TLegend> legendMassShiftEta = std::make_unique<TLegend>(0.13, 0.62, 0.7, 0.70); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
  SetLegendSettings(legendMassShiftEta.get());

  h_MmaxPTEta->SetMarkerStyle(20);
  h_MmaxPTEta->SetMarkerColor(kBlack);
  // h_MmaxPTEta->GetXaxis()->SetRangeUser(0., 17);
  h_MmaxPTEta->GetYaxis()->SetRangeUser(0.45, 0.65);
  h_MmaxPTEta->Draw("P");
  legendMassShiftEta->AddEntry(h_MmaxPTEta.get(), "#eta mass extracted", "p");
  std::unique_ptr<TLine> lineMShiftEta = std::make_unique<TLine>(0.1, 0.547862, 18, 0.547862);
  // std::unique_ptr<TLine> lineMShiftEta = std::make_unique<TLine> (1, 0.45, 3, 0.134976);
  lineMShiftEta->SetLineWidth(3);
  lineMShiftEta->SetLineColor(kAzure - 4);
  lineMShiftEta->Draw("same");
  legendMassShiftEta->AddEntry(lineMShiftEta.get(), "PDG #eta mass", "l");

  DrawingHeaderStandardLines(0.15, 0.89, 0.84, 0.79);
  lHeader->DrawLatexNDC(0.15, 0.74, "#eta #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");

  legendMassShiftEta->Draw("same");
  canvas->SaveAs("Plots_nf_data/h_MmaxPTEta.png");
  // canvas->SaveAs("/Users/lauragans-bartl/Downloads/h_MmaxPTEta.pdf");

  //**Pion
  std::unique_ptr<TLegend> legendMassShiftPion = std::make_unique<TLegend>(0.13, 0.49, 0.7, 0.57); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
  SetLegendSettings(legendMassShiftPion.get());
  h_MmaxPTPion->Draw("P");
  h_MmaxPTPion->GetYaxis()->SetRangeUser(0.123, 0.137);
  legendMassShiftPion->AddEntry(h_MmaxPTPion.get(), "#pi^{0} mass extracted", "p");
  std::unique_ptr<TLine> lineMShiftPion = std::make_unique<TLine>(0.1, 0.134976, 18, 0.134976);
  lineMShiftPion->SetLineWidth(3);
  lineMShiftPion->SetLineColor(kAzure - 4);
  lineMShiftPion->Draw("same");
  legendMassShiftPion->AddEntry(lineMShiftPion.get(), "PDG #pi^{0} mass", "l");

  DrawingHeaderStandardLines(0.15, 0.75, 0.70, 0.65);
  lHeader->DrawLatexNDC(0.15, 0.60, "#pi^{0} #rightarrow e^{+} e^{-} #gamma, rec. PCMDalitz");

  legendMassShiftPion->Draw("same");
  canvas->SaveAs("Plots_nf_data/h_MmaxPTPion.png");
  // canvas->SaveAs("/Users/lauragans-bartl/Downloads/h_MmaxPTPion.pdf");
}
