#include "helper.h"
#include "TGraph.h"
void plottingNfData()
{

  ///****PRE-SETTINGS *****/
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  // canvas
  std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
  setCanvasStandardSettings(canvas.get());

  // header
  std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>();
  setHeaderSettings(lHeader.get());

  // TFile* usedFile;

  // open file
  //  TFile* readfile = SafelyOpenRootfile("/Users/lauragans-bartl/master/data/AnalysisResults_nf_dd.root");
  // TFile* readfile = safelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC25.root");
  TFile* readfile = safelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24.root");

  std::unique_ptr<THnSparseD> hmptSame(readfile->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> hmptMix(readfile->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision(readfile->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));

  TFile* readfileTrue = safelyOpenRootfile("/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24f4dT.root");
  std::unique_ptr<THnSparseD> hmptTrue(readfileTrue->Get<THnSparseD>("pi0eta-to-gammagamma-mc-pcmdalitzee/Pair/Eta/hs_Primary"));
  std::unique_ptr<TH1F> hmptGen(readfileTrue->Get<TH1F>("pi0eta-to-gammagamma-mc-pcmdalitzee/Generated/Eta/hPt"));
  if(hmptGen==nullptr){
    std::cout << "I am empty \n";
  }

  //fuer die Funtion drawingHeaderStandardLines
  const char* usedDataset = "LHC24, LHC24f4d"; 
  const char* collisonAndBfield = "pp, #sqrt{#it{s}} = 13.6 TeV";  
  const char* particleDecayEta = "#eta #rightarrow e^{+} e^{-} #gamma";
  const char* particleDecayPion = "#pi^{0} #rightarrow e^{+} e^{-} #gamma";

  //"process" switches
  bool calculatingSOverB = false;
  bool investigatingMassShift = false;
  bool drawingHistos = true;

  //Defining Histos
  std::unique_ptr<TH1D> hEtaPt= std::make_unique<TH1D>("hEtaPt"," ; #it{p}_{T}; #it{N}",14, 0., 14.);
  setHistoStandardSettings1D(hEtaPt.get());

  // Defining vectors for projection
  std::vector<TH1D*> vechProjectSame;
  std::vector<TH1D*> vechProjectMix;
  std::vector<TH1D*> vechProjectTue;
  std::vector<std::string> vecProjType = {"minv", "pt"};

  std::vector<double> vecPtMin = {0.1, 0.4, 0.5, 0.7, 0.8, 0.9, 1., 2., 3., 5., 18.};

  // Fittting
  std::unique_ptr<TF1> gaus = std::make_unique<TF1>("gaus", "gaus", 0.47, 0.57);                           // original: 0.40, 0.60
  std::unique_ptr<TF1> gaus2 = std::make_unique<TF1>("gaus2", "gaus", 0.05, 0.17);                         // original: 0.40, 0.60
  std::unique_ptr<TF1> fSigAndBack = std::make_unique<TF1>("fSigAndBack", "gaus(0)+ pol2(3)", 0.40, 0.63); // gaus: [0]=A, [1]=mean, [2]=sigma //Zahlen in den Klammern: sagen wo die Paramter fuer die funktionen zu finden sind
  std::unique_ptr<TF1> pol2 = std::make_unique<TF1>("pol2", "pol2", 0., 0.8);

  // Massepeak in Pt
  std::unique_ptr<TH1D> hMmaxPTEta = std::make_unique<TH1D>(" hMmaxPTEta", ";#it{p}_{T} (GeV/#it{c}); #LT #it{m}_{#eta} #GT (GeV/#it{c}^{2})", vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hMmaxPTEta.get(), 1.2, 1.6);

  std::unique_ptr<TH1D> hMmaxPTPion = std::make_unique<TH1D>(" hMmaxPTPion", ";#it{p}_{T} (GeV/#it{c}); #LT #it{m}_{#pi_{0}} #GT (GeV/#it{c}^{2})", vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hMmaxPTPion.get(), 1.2, 1.7);

  //Effizienz
  std::unique_ptr<TH1D> hNEtaPt= std::make_unique<TH1D>("hNEtaPt"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaPt.get());

  std::unique_ptr<TH1D> hNEtaPtRelErr= std::make_unique<TH1D>("hNEtaPtRelErr"," ; #it{p}_{T} (GeV/#it{c}); relativ statistical uncertainty",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaPtRelErr.get());

  std::unique_ptr<TH1D> hNEtaPtMCRelErr= std::make_unique<TH1D>("hNEtaPtMCRelErr"," ; #it{p}_{T} (GeV/#it{c}); relativ statistical uncertainty",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaPtMCRelErr.get());

  std::unique_ptr<TH1D> hNEtaPtAll= std::make_unique<TH1D>("hNEtaPtAll"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaPtAll.get());

  std::unique_ptr<TH1D> hNEtaPtAllTrue= std::make_unique<TH1D>("hNEtaPtAllTrue"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaPtAllTrue.get(), 1.2, 1.4);

  std::unique_ptr<TH1D> hNEtaCompair= std::make_unique<TH1D>("hNEtaCompair"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hNEtaCompair.get());
  //OutputFile
  std::ofstream out("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/output.txt");

  // Projections aus THNSparse
  for (size_t proj = 0; proj < 2; proj++) { // 0 = minv - histo, 1 = pt - histo

    TH1D* hSame = hmptSame->Projection(static_cast<Int_t>(proj));
    hSame->SetName(Form("h_name_same_%zu", proj));
    setHistoStandardSettings1D(hSame);
    hSame->Draw("pe");
    vechProjectSame.push_back(hSame);
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_same_%s.png", vecProjType[proj].data()));

    TH1D* hMix = hmptMix->Projection(static_cast<Int_t>(proj));
    hMix->SetName(Form("h_name_mix_%zu", proj));
    setHistoStandardSettings1D(hMix);
    hMix->Draw("pe");
    vechProjectMix.push_back(hMix);
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_mix_%s.png", vecProjType[proj].data()));

    TH1D* hTrue = hmptTrue->Projection(static_cast<Int_t>(proj));
    hTrue->SetName(Form("h_name_true_%zu", proj));
    setHistoStandardSettings1D(hTrue);
    hTrue->Draw("pe");
    vechProjectTue.push_back(hTrue);
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_true_%s.png", vecProjType[proj].data()));
  }


  //****PT-Ptojections****//
  // Pt-Bins werden gesetzt
  std::vector<std::string> vecPtMinLegend = {"0.1", "0.4", "0.5", "0.7", "0.8", "0.9", "1", "2", "3", "5", "18"};
  std::vector<int> vecRebin = {8, 8, 5, 5, 5, 5, 5, 5, 8, 8};
  double epsilon = 1.e-9; // Hilfsmittel um sicherzustellen, dass der richtige Bin gelesen wird
  std::vector<int> vecPtMinBinSame;
  std::vector<int> vecPtMinBinMix;
  std::vector<int> vecPtMinBinTrue;
  // Passende Bin Nummern werden aus dem pt Histogramm geholt
  for(const auto& ptMin : vecPtMin){      //TO DO: WIESO DREI MAL???
    int iSame = vechProjectSame[1]->FindBin(ptMin + epsilon); // epsilon damit die richtigen Bins gelesen werden
    vecPtMinBinSame.push_back(iSame);

    int iMix = vechProjectMix[1]->FindBin(ptMin + epsilon);
    vecPtMinBinMix.push_back(iMix);

    int iTrue = vechProjectTue[1]->FindBin(ptMin + epsilon);
    vecPtMinBinTrue.push_back(iTrue);
  }

  // int nRebin = 5;
  // projection der Masseverteilung in den einzelnen pt Bereichen
  std::vector<TH1D*> vechSameProjmPtranges;
  std::vector<TH1D*> vechMixProjmPtranges;
  std::vector<TH1D*> vechTrueProjmPtranges;
  std::vector<double> vecNGenEtaPerPt;
  std::vector<double> vecNGenEtaPerPtError;
  for (size_t i = 0; i < vecPtMinBinSame.size() - 1; i++) { // Projections of specific pt-ranges
    // same
    auto* hmptSameClone = (THnSparseD*)hmptSame->Clone(Form("hmpt_same_clone_%zu", i));
    hmptSameClone->GetAxis(1)->SetRange(vecPtMinBinSame[i], vecPtMinBinSame[i + 1]); // SetRange erwartet Bins. Falls keine Bins SetRangeUsers
    TH1D* hSame = hmptSameClone->Projection(0);                                             // masseverteilung wird rausgeholt
    hSame->SetName(Form("h_name_same_projection_%zu", i));
    setHistoStandardSettings1D(hSame);
    hSame->Rebin(vecRebin[i]); // kombiniert einzelne Bins zusammen. so viele wie in der Klammer steht.
    hSame->Draw("pe");
    vechSameProjmPtranges.push_back(hSame);
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_same_%0.1f-%0.1f_Rebin%i.png", vecPtMin[i], vecPtMin[i+1], nRebin));

    // mixed
    auto* hmptMixClone = (THnSparseD*)hmptMix->Clone(Form("hmptMixClone%zu", i));
    hmptMixClone->GetAxis(1)->SetRange(vecPtMinBinMix[i], vecPtMinBinMix[i + 1]);
    TH1D* hMix = hmptMixClone->Projection(0);
    hMix->SetName(Form("h_name_mix_projection_%zu", i));
    setHistoStandardSettings1D(hMix);
    hMix->Rebin(vecRebin[i]);
    hMix->Draw("pe");
    vechMixProjmPtranges.push_back(hMix);
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_mix_%0.1f-%0.1f_Rebin%i.png", vecPtMin[i], vecPtMin[i+1], nRebin));

    // true
    auto* hmptTrueClone = (THnSparseD*)hmptTrue->Clone(Form("hmptTrueClone%zu", i));
    hmptTrueClone->GetAxis(1)->SetRange(vecPtMinBinTrue[i], vecPtMinBinTrue[i + 1]);
    TH1D* hTrue = hmptTrueClone->Projection(0);
    hTrue->SetName(Form("h_name_true_projection_%zu", i));
    setHistoStandardSettings1D(hTrue);
    hTrue->Rebin(vecRebin[i]);
    hTrue->Draw("pe");
    vechTrueProjmPtranges.push_back(hTrue);
    std::cout << "vecPtMinBinTrue[i]" << vecPtMinBinTrue[i] << " | vecPtMinBinTrue[i + 1]:" << vecPtMinBinTrue[i + 1] << "\n";
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hproj_true_%0.1f-%0.1f_Rebin%i.png", vecPtMin[i], vecPtMin[i+1], nRebin));
    double errTrue = 0;
    double nrecTrueEtas = hTrue->IntegralAndError(hTrue->FindBin(0.45), hTrue->FindBin(0.6),errTrue);
    hNEtaPtMCRelErr->SetBinContent(static_cast<Int_t>(i) + 1, 100 * errTrue/nrecTrueEtas);

    //MC: true Eta / generated
    double err = 0;
    double nGenEtas = hmptGen->IntegralAndError(vecPtMinBinTrue[i], vecPtMinBinTrue[i + 1],err);
    vecNGenEtaPerPt.push_back(nGenEtas);
    vecNGenEtaPerPtError.push_back(err);
    
  }

  // // Getting Binning of original histos:
  // std::vector<double> vecPtBins;
  // auto nBins = static_cast<size_t>(vechSameProjmPtranges[0]->GetNbinsX());
  // vecPtBins.resize(nBins + 1);
  // for (size_t iBin = 1; iBin <=  nBins + 1; iBin++) {
  //   const int bin = static_cast<int>(iBin) + 1; 
  //   vecPtBins.at(iBin - 1) = vechSameProjmPtranges[0]->GetBinLowEdge(bin);
  // }

  //****Signal Extraction****//
  // Defining Vectors for next loop
  std::vector<double> vecEtaSoverBMixed;
  std::vector<double> vecPionSoverBMixed;

  for (size_t i = 0; i < vechSameProjmPtranges.size(); i++) {
    //**1.1 Defining Histos
    TH1D* hSame = (TH1D*)vechSameProjmPtranges[i]->Clone(Form("h_same_%zu", i));
    TH1D* hMix = (TH1D*)vechMixProjmPtranges[i]->Clone(Form("h_mix_%zu", i));
    TH1D* hTrue = (TH1D*)vechTrueProjmPtranges[i]->Clone(Form("h_true_%zu", i));
    //**1.2 Scaling
    //**1.2.1 Scaling Version 1 - one scaling factor
    // double skalar_same = hSame->Integral(hSame->FindBin(0.7), hSame->FindBin(0.8));
    // double skalar_mix = hMix->Integral(hMix->FindBin(0.7), hMix->FindBin(0.8));
    // double scaling_factor = skalar_same / skalar_mix;

    //**1.2.2 scaling Version 2 - scaling function
    TH1D* hRatioSameMix = makingRatio(hSame, hMix, "hRatioSameMix");
    int binMin = hRatioSameMix->FindBin(0.44);
    int binMax = hRatioSameMix->FindBin(0.58);
    for (int iBin = binMin; iBin <= binMax; iBin++) {
      hRatioSameMix->SetBinContent(iBin, 0);
      hRatioSameMix->SetBinError(iBin, 0);
    }
    hRatioSameMix->SetAxisRange(0.3, 0.8);
    hRatioSameMix->GetYaxis()->SetTitle("same/mix");
    hRatioSameMix->Draw("pe");
    hRatioSameMix->Fit(pol2.get(), "0M", "", 0.3, 0.8);
    //**1.3 Drawing Scaled mixed Histo
    TH1D* hMixScaled = (TH1D*)vechMixProjmPtranges[i]->Clone(Form("h_mix_scaled_%zu", i));
    setHistoStandardSettings1D(hMixScaled);
    // hMixScaled->Scale(scaling_factor);         //scaling factor
    hMixScaled->Multiply(pol2.get()); // scaling function
    hMixScaled->Draw("pe");
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/h_mix_scaled_%0.1f-%0.1f.png",vecPtMin[i], vecPtMin[i+1]));

    TH1* hSameScaledMixRatio = makingRatio(hSame, hMixScaled, "hSameScaledMixRatio");
    hSameScaledMixRatio->Draw("");
    canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hSameScaledMixRatio%0.1f-%0.1f.png", vecPtMin[i], vecPtMin[i + 1]));

    // TH1* hSameMixRatio = makingRatio(hSame, hMix, "hSameMixRatio");
    // hSameMixRatio->Draw("");
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hSameMixRatio_%0.1f-%0.1f.png",vecPtMin[i], vecPtMin[i+1]));

    //** 2. Calculating Eta-Signal
    std::unique_ptr<TLegend> legendMass = std::make_unique<TLegend>(0.13, 0.49, 0.35, 0.61); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendMass.get());
    TH1D* hSignal = (TH1D*)vechSameProjmPtranges[i]->Clone(Form("h_same_sig_%zu", i));
    hSignal->Add(hSignal, hMixScaled, 1., -1.);
    hSignal->GetYaxis()->SetTitle("#it{N}");
    hSignal->GetXaxis()->SetTitle("#it{m}_{ee#gamma} (GeV/#it{c}^{2})");
    setHistoStandardSettings1D(hSignal);

    hSignal->Fit(gaus.get(), "0M", "", 0.470, 0.570);
    gaus->SetLineColor(kRed);
    //**Fitting Eta
    fSigAndBack->SetParLimits(1, 0.95 * gaus->GetParameter(1), 1.05 * gaus->GetParameter(1));
    fSigAndBack->SetParLimits(2, 0.5 * gaus->GetParameter(2), 1.05 * gaus->GetParameter(2));
    fSigAndBack->SetLineColor(kRed + 1);
    fSigAndBack->SetLineWidth(3);
    fSigAndBack->SetNpx(10000);
    hSignal->Fit(fSigAndBack.get(), "0M", "", 0.2, 0.8);
    // hSignal->GetXaxis()->SetRangeUser(0.28, 0.7);

    //**Einschub: Massepeaks werden rausgezogen, em die Abweichung von der erwareten Masse zu plotten. Plot: nach dem Loop
    double mMaxEta = gaus->GetParameter(1);
    double mMaxEtaError = gaus->GetParError(1);
    hMmaxPTEta->SetBinContent(static_cast<Int_t>(i) + 1, mMaxEta);
    hMmaxPTEta->SetBinError(static_cast<Int_t>(i) + 1, mMaxEtaError);

    //**Fitting Pion with Steffis function
    TF1* fitPion = fitPi0ShapeSave(hSignal, 0.02, 0.25);
    fitPion->SetLineColor(kBlue);
    fitPion->SetLineWidth(3);
    fitPion->Draw("same");

    //**Einschub: Massepeaks werden rausgezogen, em die Abweichung von der erwareten Masse zu plotten. Plot: nach dem Loop
    double mMaxPion = fitPion->GetParameter(4);
    double mMaxPionError = fitPion->GetParError(4);

    hMmaxPTPion->SetBinContent(static_cast<Int_t>(i) + 1, mMaxPion);
    hMmaxPTPion->SetBinError(static_cast<Int_t>(i) + 1, mMaxPionError);

    hSignal->Draw("pe");
    hSignal->GetYaxis()->SetRangeUser(-300., 6000);
    fSigAndBack->Draw("same");

    //Berechnen der Effizienz
    // hNCollision->Draw();
    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNCollision.png"));
    double_t nCollisions = hNCollision->GetBinContent(12);

    double mean = fSigAndBack->GetParameter(1);
    double sigma = fSigAndBack->GetParameter(2);

    double xmin = mean - (3.0*sigma);
    double xmax = mean + (3.0*sigma);

    int bmin = hSignal->GetXaxis()->FindBin(xmin);
    int bmax = hSignal->GetXaxis()->FindBin(xmax);
    double err = 0;

    double nEtas = hSignal->IntegralAndError(bmin, bmax, err);
    double nEtasPColl = nEtas/nCollisions;
    double errPColl = err/nCollisions;

    double nTrueEtas = hTrue->GetEntries();
    // double nTruEteasPColl = nTrueEtas/nCollisions;

    
    hNEtaPt->SetBinContent(static_cast<Int_t>(i) + 1, nEtasPColl);
    hNEtaPt->SetBinError(static_cast<Int_t>(i) + 1, errPColl);

    hNEtaPtRelErr->SetBinContent(static_cast<Int_t>(i) + 1, 100 * (errPColl/nEtasPColl));
    // hNEtaPtMCRelErr->SetBinContent(static_cast<Int_t>(i) + 1, 100 * (errPColl/nEtasPColl));
    
    double nEtaoBR = nEtas/(nCollisions * 0.007);
    hNEtaPtAll->SetBinContent(static_cast<Int_t>(i) + 1, nEtaoBR);
    
    double nTrueEtaoBR = nTrueEtas/(nCollisions * 0.007);
    hNEtaPtAllTrue->SetBinContent(static_cast<Int_t>(i) + 1, nTrueEtaoBR);
    
    out << "\n" <<nCollisions << "\n" << "Pt-Slice " << i << ": " << nEtas << " and nEtasPColl:" << nEtasPColl << " (with Mean:" << mean << " and sigma:" << sigma << " )" << "\n" << "Trues: " << nTrueEtas << "with: " << nTrueEtaoBR <<  "\n";

    //MC: true Eta / generated
    double nGenEta = vecNGenEtaPerPt[i];
    if(nGenEta==0){
      nGenEta=1;
    }
    double nEtasCompair = nTrueEtas/nGenEta;
    hNEtaCompair->SetBinContent(static_cast<Int_t>(i) + 1, nEtasCompair);
    
    // drawLine(0.06, -500, 0.06, 2000, 3, 800);
    // drawLine(0.18, -500, 0.18, 2000, 3, 800);
    // drawLine(0.46, -500, 0.46, 2000, 3, 800);
    // drawLine(0.58, -500, 0.58, 2000, 3, 800);


    // hSignal->GetYaxis()->SetRangeUser(-100., 50e3);
    drawingMarker(hSame, kGray);
    hSame->SetMarkerColor(kGray);
    hSame->SetLineColor(kGray);
    // hSame->Draw("same");

    hTrue->Draw("same");
    // hTrue->Scale(3.3);
    drawingMarker(hTrue, 47, 2, kOrange+8);

    legendMass->AddEntry(hSignal, "signal, data", "p");
    // legendMass->AddEntry(hSame, "same event", "p");
    legendMass->AddEntry(hTrue, "true #eta, MC", "p");
    // legendMass->AddEntry(fitPion, "Crystall ball fit", "l");
    legendMass->AddEntry(fSigAndBack.get(), "pol2 + Gausian fit", "l");

    legendMass->Draw("same");
    canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/h_signal_eta_pion%0.1f-%0.1f.png", vecPtMin[i], vecPtMin[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/h_signal_eta_pion%0.1f-%0.1f.pdf", vecPtMin[i], vecPtMin[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/h_signal_eta_pion_%0.1f-%0.1f.pdf", vecPtMin[i], vecPtMin[i+1]));        //fuer den Mas

    // hTrue->GetYaxis()->SetTitle("#it{N}");
    // hTrue->GetYaxis()->SetTitleOffset(static_cast<float>(1.6));
    // hTrue->GetXaxis()->SetTitle("#it{m}_{ee#gamma}");
    drawingMarker(hTrue);
    // hTrue->Draw("pe");

    // drawingHeaderStandardLines(0.13, usedDataset, collisonAndBfield, particleDecayEta);
    // lHeader->DrawLatexNDC(0.13, 0.65, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vecPtMinLegend[i].data(), vecPtMinLegend[i+1].data()));

    // canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hTrue%0.1f-%0.1f.png", vecPtMin[i], vecPtMin[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/hTrue%0.1f-%0.1f.pdf", vecPtMin[i], vecPtMin[i+1]));        //fuer den Mas



    //**NUR Pion Signal
    hSignal->GetXaxis()->SetRangeUser(0., 0.2);

    canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hSignalPion_%0.1f_%0.1f.png", vecPtMin[i], vecPtMin[i + 1]));
    // canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/hSignalPion_%0.1f-%0.1f.pdf", vecPtMin[i], vecPtMin[i+1]));

    hSignal->GetXaxis()->SetRangeUser(0.3, 0.7);
    hSignal->GetYaxis()->SetRangeUser(-500, 2000); // fuer den kleinen Eta Bereich

    //**NUR ETA Signal
    // Legende fuer 0.8< pt < 0.9
    drawingHeaderStandardLines(0.13, 0.85, usedDataset, collisonAndBfield, particleDecayEta);
    lHeader->DrawLatexNDC(0.13, 0.65, Form("%s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}", vecPtMinLegend[i].data(), vecPtMinLegend[i + 1].data()));

    canvas->SaveAs(Form("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hSignalEta_%0.1f_%0.1f.png", vecPtMin[i], vecPtMin[i + 1]));
    canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/Plots/NfData/hSignalEta_%0.1f-%0.1f.pdf", vecPtMin[i], vecPtMin[i+1]));

    //**3.Signal/Background (Wert fuer den jeweiligen pt-Bereich)
    //**3.1: fuer Eta
    double etaSignal = hSignal->Integral(hSignal->FindBin(0.5), hSignal->FindBin(0.55)); //->FindBin(0.46), hSignal->FindBin(0.57));

    //**Version1: scaled mixed als Hintergrund
    double etaBackgroundMixed = hMixScaled->Integral(hMixScaled->FindBin(0.5), hMixScaled->FindBin(0.55));
    double sOverBMixedEta = etaSignal / etaBackgroundMixed;
    vecEtaSoverBMixed.push_back(sOverBMixedEta);

    //**Version2: same Signal als Hintergrund
    // hSame->Add(hSame, hSignal, 1, -1);
    // double EtaBackgroundSameSignal = hSame->Integral(hSame->FindBin(0.5), hSame->FindBin(0.55));
    // double EtaSoverB_same_signal = etaSignal/EtaBackgroundSameSignal;

    //**3.2 fuer Pion
    double pionSignal = hSignal->Integral(hSignal->FindBin(0.08), hSignal->FindBin(0.14));

    //**Version1: scaled mixed als Huntergrund
    double pionBackgroundMixed = hMixScaled->Integral(hMixScaled->FindBin(0.08), hMixScaled->FindBin(0.14));
    double pionSoverBMixed = pionSignal / pionBackgroundMixed;
    vecPionSoverBMixed.push_back(pionSoverBMixed);

    //**Version2: same-etaSignal als Hintergrund
    // hSame->Add(hSame, hSignal, 1, -1);
    // double pionBackgroundSameSignal = hSame->Integral(hSame->FindBin(0.6), hSame->FindBin(0.18));
    // double PionSoverB_same_signal = pionSignal/pionBackgroundSameSignal;
  }

  if(calculatingSOverB){
    //**4.4 Signal/Background ****//
    //**4.4.1 Eta
    canvas->SetLogy(1);
    std::unique_ptr<TH1D> hEtaSoverB = std::make_unique<TH1D>("hEtaSoverB", " ;#it{p}_{T} (GeV/#it{c}); S/B", vecPtMin.size() - 1, vecPtMin.data());
    setHistoStandardSettings1D(hEtaSoverB.get());
    hEtaSoverB->GetXaxis()->SetRangeUser(0.8, 6);
    hEtaSoverB->GetYaxis()->SetRangeUser(1e-2, 12);
    for (size_t i = 0; i < vecPtMin.size(); i++) {
      if (vecPtMin[i] > 0.8) {
        hEtaSoverB->SetBinContent(hEtaSoverB->GetXaxis()->FindBin(vecPtMin[i]), vecEtaSoverBMixed[i]);
        hEtaSoverB->SetBinError(hEtaSoverB->GetXaxis()->FindBin(vecPtMin[i]), 10e-10);
      }
    }
    hEtaSoverB->SetMarkerColor(kRed + 1);
    hEtaSoverB->SetLineColor(kRed + 1);

    hEtaSoverB->Draw("");

    // canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hEtaSoverB.png");
    // canvas->SaveAs("/Users/lauragans-bartl/Downloads/hEtaSoverB.pdf");

    //**4.4.2 Pion
    std::unique_ptr<TH1D> hPionSoverB = std::make_unique<TH1D>("hPionSoverB", " ;#it{p}_{T} (GeV/#it{c}); S/B", vecPtMin.size() - 1, vecPtMin.data());
    setHistoStandardSettings1D(hPionSoverB.get());
    for (size_t i = 0; i < vecPtMin.size() - 1; i++) {
      hPionSoverB->SetBinContent(static_cast<Int_t>(i) + 1, vecPionSoverBMixed[i]);
      hPionSoverB->SetBinError(static_cast<Int_t>(i) + 1, 10e-10);
    }
    hPionSoverB->SetMarkerColor(kBlue);
    hPionSoverB->SetLineColor(kBlue);
    hPionSoverB->Draw("same");

    // canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hPionSoverB.png");

    //**Beide auf einem Canvas
    std::unique_ptr<TLegend> legendSoverB = std::make_unique<TLegend>(0.36, 0.24, 0.65, 0.32); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendSoverB.get());
    legendSoverB->AddEntry(hPionSoverB.get(), "#pi^{0}: 0.06 GeV/#it{c}^{2} < #it{m}_{ee#gamma} < 0.18 GeV/#it{c}^{2}", "p");
    legendSoverB->AddEntry(hEtaSoverB.get(), "#eta : 0.5 GeV/#it{c}^{2} < #it{m}_{ee#gamma} < 0.55 GeV/#it{c}^{2}", "p");
    
    drawingHeaderStandardLines(0.4, 0.50, usedDataset, collisonAndBfield, particleDecayEta);

    // h_PionSoverB_signal->Draw("same");
    legendSoverB->Draw("same");
    // canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/h_PionEtaSoverB.png");

    canvas->SetLogy(0);
  }

  if(investigatingMassShift){
    //***Massverschiebung
    //**Eta
    canvas->SetLeftMargin(static_cast<float>(0.11));
    std::unique_ptr<TLegend> legendMassShiftEta = std::make_unique<TLegend>(0.13, 0.62, 0.7, 0.70); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendMassShiftEta.get());

    hMmaxPTEta->SetMarkerStyle(20);
    hMmaxPTEta->SetMarkerColor(kBlack);
    // hMmaxPTEta->GetXaxis()->SetRangeUser(0., 17);
    hMmaxPTEta->GetYaxis()->SetRangeUser(0.45, 0.65);
    hMmaxPTEta->Draw("P");
    legendMassShiftEta->AddEntry(hMmaxPTEta.get(), "#eta mass extracted", "p");
    std::unique_ptr<TLine> lineMShiftEta = std::make_unique<TLine>(0.1, 0.547862, 18, 0.547862);
    // std::unique_ptr<TLine> lineMShiftEta = std::make_unique<TLine> (1, 0.45, 3, 0.134976);
    lineMShiftEta->SetLineWidth(3);
    lineMShiftEta->SetLineColor(kAzure - 4);
    lineMShiftEta->Draw("same");
    legendMassShiftEta->AddEntry(lineMShiftEta.get(), "PDG #eta mass", "l");

    drawingHeaderStandardLines(0.15, 0.89, usedDataset, collisonAndBfield, particleDecayEta);

    legendMassShiftEta->Draw("same");
    // canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hMmaxPTEta.png");
    // canvas->SaveAs("/Users/lauragans-bartl/Downloads/hMmaxPTEta.pdf");

    //**Pion
    std::unique_ptr<TLegend> legendMassShiftPion = std::make_unique<TLegend>(0.13, 0.49, 0.7, 0.57); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendMassShiftPion.get());
    hMmaxPTPion->Draw("P");
    hMmaxPTPion->GetYaxis()->SetRangeUser(0.123, 0.137);
    legendMassShiftPion->AddEntry(hMmaxPTPion.get(), "#pi^{0} mass extracted", "p");
    std::unique_ptr<TLine> lineMShiftPion = std::make_unique<TLine>(0.1, 0.134976, 18, 0.134976);
    lineMShiftPion->SetLineWidth(3);
    lineMShiftPion->SetLineColor(kAzure - 4);
    lineMShiftPion->Draw("same");
    legendMassShiftPion->AddEntry(lineMShiftPion.get(), "PDG #pi^{0} mass", "l");

    drawingHeaderStandardLines(0.15, 0.75, usedDataset, collisonAndBfield, particleDecayPion);

    legendMassShiftPion->Draw("same");
    // canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hMmaxPTPion.png");
    // canvas->SaveAs("/Users/lauragans-bartl/Downloads/hMmaxPTPion.pdf");
  }
  
  //****Drawing****//
  if(drawingHistos){
    //** hNEtaPt
    hNEtaPt->GetYaxis()->SetTitle("NEta/collision");
    hNEtaPt->Draw("P");
    hNEtaPt->SetMarkerStyle(20);
    canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNEtaPt.png");

    //** hNEtaPtRelErr und legendEtaEffiErr
    std::unique_ptr<TLegend> legendEtaEffiErr = std::make_unique<TLegend>(0.38, 0.70, 0.7, 0.8); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendEtaEffiErr.get());
    hNEtaPtRelErr->GetYaxis()->SetTitle("relative statistical uncertainty");
    hNEtaPtRelErr->Draw("P");
    hNEtaPtRelErr->SetMarkerStyle(20);
    hNEtaPtRelErr->GetXaxis()->SetRangeUser(1, 5);
    hNEtaPtRelErr->GetYaxis()->SetRangeUser(0, 15);
    legendEtaEffiErr->AddEntry(hNEtaPtRelErr.get(), "data", "p");
    legendEtaEffiErr->AddEntry(hNEtaPtMCRelErr.get(), "true MC", "p");
    hNEtaPtMCRelErr->SetMarkerStyle(24);
    hNEtaPtMCRelErr->Draw("same, p");
    legendEtaEffiErr->Draw();
    canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNEtaPtRelErr.png");
    canvas->SaveAs("/Users/lauragans-bartl/Downloads/Plots/NfData/hNEtaPtRelErr.png");
    std::cout << "Printing hNEtaPtMCRelErr\n";
    TGraph gr(hNEtaPtMCRelErr.get());
    gr.Print();


    //** hNEtaPtAll
    std::unique_ptr<TLegend> legendEtaEffi = std::make_unique<TLegend>(0.38, 0.76, 0.7, 0.8); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendEtaEffi.get());
    hNEtaPtAll->GetYaxis()->SetTitle("NEta/(collision * BR)");
    hNEtaPtAll->Draw("P");
    hNEtaPtAll->GetXaxis()->SetRangeUser(0, 5.);
    hNEtaPtAll->GetYaxis()->SetRangeUser(-1e-6, 50e-6);
    hNEtaPtAll->SetMarkerStyle(20);
    legendEtaEffi->AddEntry(hNEtaPtAll.get(), "#eta meson", "p");
    lHeader->DrawLatexNDC(0.4, 0.85, "ALICE work in progress");
    lHeader->DrawLatexNDC(0.4, 0.80, "LHC25ah");
    legendEtaEffi->Draw("same");
    canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNEtaPtAll.png"); 
    canvas->SaveAs("/Users/lauragans-bartl/Downloads/hNEtaPtAll.pdf");
    
    //** hNEtaPtAllTrue
    std::unique_ptr<TLegend> legendEtaEffiTrue = std::make_unique<TLegend>(0.38, 0.76, 0.7, 0.8); // RELATIV!!!(x1,y1,x2,y2), HINTER den Draw Befehl von der FUnktion, weil die sonst alles übermalt
    setLegendSettings(legendEtaEffiTrue.get());
    hNEtaPtAllTrue->Draw("P");
    hNEtaPtAllTrue->GetXaxis()->SetRangeUser(0, 5.);
    hNEtaPtAllTrue->GetYaxis()->SetTitle("NTrueEta/(collision * BR)");
    hNEtaPtAllTrue->SetMarkerStyle(20);
    legendEtaEffiTrue->AddEntry(hNEtaPtAllTrue.get(), "#eta meson", "p");
    
    lHeader->DrawLatexNDC(0.4, 0.85, "ALICE work in progress");
    lHeader->DrawLatexNDC(0.4, 0.80, "LHC25ah");
    legendEtaEffiTrue->Draw();
    canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNEtaPtAllTrue.png");
    canvas->SaveAs("/Users/lauragans-bartl/Downloads/hNEtaPtAllTrue.pdf");

    //** hNEtaCompair
    hNEtaCompair->GetYaxis()->SetTitle("True Eta / generated Eta");
    hNEtaCompair->Draw("P");
    hNEtaCompair->SetMarkerStyle(20);
    canvas->SaveAs("/Users/lauragans-bartl/master/MyAnalysis/Plotting/Plots/Plots_nf_data/hNEtaCompair.png");
  }

  out.close();

}
int main(){
  plottingNfData();
}
