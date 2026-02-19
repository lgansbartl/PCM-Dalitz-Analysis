///////////////////////////////////////////////////////////////////////////////////
// Macro for SIGNAL EXTRACTION 
///////////////////////////////////////////////////////////////////////////////////

//To do:
// true Information mit einbauen -> wie umgehen, wenn es keine true gibt?

//******************INCLUDINGS******************//
#include "helperPlotting.h"


void extractSignal(const char* filename, 
                  const char* filenameMC = "",
                  const std::string filenameSaving = "savingExtract.root", 
                  bool isMC = false, 
                  double minv = 0.55,
                  bool debug = true)                 
{

  ////*****FILES****/////
  std::ofstream out("debugFile.txt");

  TFile* fileRecMC = safelyOpenRootfile(filename);
  std::unique_ptr<THnSparseD> sSame(fileRecMC->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/same/hs"));
  std::unique_ptr<THnSparseD> sMix(fileRecMC->Get<THnSparseD>("pi0eta-to-gammagamma-pcmdalitzee/Pair/mix/hs"));
  std::unique_ptr<TH1D> hNCollision(fileRecMC->Get<TH1D>("pi0eta-to-gammagamma-pcmdalitzee/Event/after/hCollisionCounter"));
  hNCollision->SetName("hNCollision");
  double_t nCollisions = hNCollision->GetBinContent(12);
  

  
  std::unique_ptr<THnSparseD> sTrue;
  std::unique_ptr<TH1F> hGen;

  TFile* fileTrueMC = safelyOpenRootfile(filenameMC);

  std::vector<double> vecOriginRanges;
  if(isMC){
    auto* obj = fileTrueMC->Get<THnSparseD>("pi0eta-to-gammagamma-mc-pcmdalitzee/Pair/Eta/hs_Primary");
    sTrue.reset(static_cast<THnSparseD*>(obj->Clone("hsPrimaryClone")));

    auto* obj2 = fileTrueMC->Get<TH1F>("pi0eta-to-gammagamma-mc-pcmdalitzee/Generated/Eta/hPt");
    hGen.reset(static_cast<TH1F*>(obj2->Clone("hGen")));
    setHistoStandardSettings1D(hGen.get());
    

    if(debug){
      // TH1D* hSame = sSame->Projection(static_cast<Int_t>(1));
      for (int i = 1; i <= hGen->GetNbinsX(); ++i) {
        double x1 = hGen->GetXaxis()->GetBinLowEdge(i);
        double x2 = hGen->GetXaxis()->GetBinUpEdge(i);
        out << "x1: " << x1 << "| x2: " << x2 << "\n";
        vecOriginRanges.push_back(x1);
      }
    }
  }
  
  
  std::unique_ptr<TFile> savingFile = std::unique_ptr<TFile>(TFile::Open(filenameSaving.c_str(), "RECREATE"));    //.c_str() casts std::string to char*
  hNCollision->Write();

  ////**********SAFTY CHECKS**********/////
  if (!savingFile || !savingFile->IsOpen()) {
    std::cerr << "ERROR can't open saving file \n";
  }
  savingFile->cd();


  ////**********DEFINITIONS**********/////
  std::vector<double> vecPtMin = { 0.1, 0.2, 0.4, 0.8, 1.2, 2., 5., 18.};
  std::vector<int> vecRebin = {8, 8, 5, 5, 5, 5, 8, 8, 8};

  

  ////*****HISTOS****/////
  std::unique_ptr<TH1D> hRawpT= std::make_unique<TH1D>("hRawpT"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hRawpT.get());

  std::unique_ptr<TH1D> hTrueRawpT= std::make_unique<TH1D>("hTrueRawpT"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hTrueRawpT.get());

  std::unique_ptr<TH1D> hGenRawpT= std::make_unique<TH1D>("hGenRawpT"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hGenRawpT.get());

  std::unique_ptr<TH1D> hGenRawpTDalitz= std::make_unique<TH1D>("hGenRawpTDalitz"," ; #it{p}_{T} (GeV/#it{c}); #it{N}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hGenRawpTDalitz.get());

  std::unique_ptr<TH1D> hSigni= std::make_unique<TH1D>("hSigni"," ; #it{p}_{T} (GeV/#it{c}); #frac{S}{S+B}",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hSigni.get());

  std::unique_ptr<TH1D> hSB= std::make_unique<TH1D>("hSB"," ; #it{p}_{T} (GeV/#it{c}); S/B",vecPtMin.size() - 1, vecPtMin.data());
  setHistoStandardSettings1D(hSB.get());


  if(isMC){
    TH1D* hGenNB = (TH1D*) hGen->Rebin(static_cast<Int_t>(vecPtMin.size() - 1), "hGenNB",vecPtMin.data()) ;      
    normalizeSepc(hGenNB, nCollisions);
    hGenNB->Write();
    hGen->Write();
  }




  ////*****FITS****/////
  std::unique_ptr<TF1> fScaling = std::make_unique<TF1>("fScaling", "pol2", 0., 0.8);
  std::unique_ptr<TF1> fGaus = std::make_unique<TF1>("fGaus", "gaus", 0.47, 0.57);      
  

  ////**********PROJECTING Pt SLICES**********/////
  double epsilon = 1.e-9;
  std::vector<int> vecRangeLow;
  std::vector<int> vecRangeUp;

  ////*****Translating pt values in Bins****/////
  //******************************************************TO DO: irgendwas stimmt mit den beiden Projektionen nicht!! Fehler: "Warning in <TFile::Append>: Replacing existing TH1: Pair/same/hs_proj_1 (Potential memory leak)." */
  for(const auto& ptMin : vecPtMin){   
    int iSameUp = sSame->Projection(1)->FindBin(ptMin + epsilon);     //untere Grenze
    int iSameDown = sSame->Projection(1)->FindBin(ptMin - epsilon);     //obere Grenze
    vecRangeLow.push_back(iSameUp);
    vecRangeUp.push_back(iSameDown);

  }

  ////*****Projecting pT slides****/////
  bool rebin = true;
  std::vector<TH1D*> vecHSamePtSlides;
  std::vector<TH1D*> vecHMixPtSlides;
  std::vector<TH1D*> vecHTruePtSlides;
  std::vector<double> vecNGenEtaPerPt;
  std::vector<double> vecNGenEtaPerPtError;
  for (size_t i = 0; i < vecRangeLow.size() - 1; i++) {
    //for same
    auto* hmptSameClone = (THnSparseD*)sSame->Clone(Form("hmptSameClone%zu", i));
    hmptSameClone->GetAxis(1)->SetRange(vecRangeLow[i], vecRangeUp[i + 1]);
    TH1D* hSame = hmptSameClone->Projection(0);                                             // masseverteilung wird rausgeholt
    hSame->SetName(Form("hSameProjection%zu", i));
    setHistoStandardSettings1D(hSame);
    if(rebin) hSame->Rebin(vecRebin[i]);
    hSame->Write();
    vecHSamePtSlides.push_back(hSame);

    //for mixed
    auto* hmptMixClone = (THnSparseD*)sMix->Clone(Form("hmptMixClone%zu", i));
    hmptMixClone->GetAxis(1)->SetRange(vecRangeLow[i], vecRangeUp[i + 1]);
    TH1D* hMix = hmptMixClone->Projection(0);                                             // masseverteilung wird rausgeholt
    hMix->SetName(Form("hMixProjection%zu", i));
    setHistoStandardSettings1D(hMix);
    if(rebin) hMix->Rebin(vecRebin[i]);
    hMix->Write();
    vecHMixPtSlides.push_back(hMix);

    if(isMC) {
      auto* hmptTrueClone = (THnSparseD*)sTrue->Clone(Form("hmptTrueClone%zu", i));
      hmptTrueClone->GetAxis(1)->SetRange(vecRangeLow[i], vecRangeUp[i + 1]);
      TH1D* hTrue = hmptTrueClone->Projection(0);                                             // masseverteilung wird rausgeholt
      hTrue->SetName(Form("hTrueProjection%zu", i));
      setHistoStandardSettings1D(hTrue);
      if(rebin) hTrue->Rebin(vecRebin[i]);
      hTrue->Write();
      vecHTruePtSlides.push_back(hTrue);

      //Generated Eta
      double err = 0;
      double nGenEtas = hGen->IntegralAndError(vecRangeLow[i], vecRangeUp[i + 1],err);
      vecNGenEtaPerPt.push_back(nGenEtas);
      vecNGenEtaPerPtError.push_back(err);
    }

  }

  ////**********EXTRACTING For EACH PT SLICE**********/////
  std::vector<double> vecNEtas;
  std::vector<double> vecNEtasErr;
  for(size_t iSlide = 0; iSlide < vecHSamePtSlides.size(); iSlide++){
    TH1D* hSame = (TH1D*)vecHSamePtSlides[iSlide]->Clone(Form("hSame_%zu", iSlide));
    TH1D* hMix = (TH1D*)vecHMixPtSlides[iSlide]->Clone(Form("hMix%zu", iSlide));
    TH1D* hTrue;
    if(isMC) {
      hTrue = (TH1D*)vecHTruePtSlides[iSlide]->Clone(Form("hTrue%zu", iSlide));
    }


    ////*****SCALING****/////
    //with Scaling function
    TH1D* hRatioSameMix = makingRatio(hSame, hMix, Form("hRatioSameMix%zu", iSlide));
    double xMin = minv - 0.05;
    double xMax = minv + 0.05;
    int binMin = hRatioSameMix->FindBin(xMin);
    int binMax = hRatioSameMix->FindBin(xMax);
    for (int iBin = binMin; iBin <= binMax; iBin++) {
      hRatioSameMix->SetBinContent(iBin, 0);
      hRatioSameMix->SetBinError(iBin, 0);
    }
    hRatioSameMix->SetAxisRange(0.3, 0.8);
    hRatioSameMix->GetYaxis()->SetTitle("same/mix");
    setHistoStandardSettings1D(hRatioSameMix);
    hRatioSameMix->Draw("pe");
    hRatioSameMix->Fit(fScaling.get(), "0M", "", 0.3, 0.8);
    // hRatioSameMix->Write();


    ////*****SCALING hMIX*****/////
    TH1D* hMixScaled = (TH1D*)hMix->Clone( Form("hMixScaled%zu", iSlide));
    hMixScaled->Multiply(fScaling.get()); 
    setHistoStandardSettings1D(hMixScaled);
    hMixScaled->Write();

    TH1* hSameScaledMixRatio = makingRatio(hSame, hMixScaled, Form("hSameScaledMixRatio%zu", iSlide));
    hSameScaledMixRatio->Draw("");

    ////*****EXTRACTING SIGNAL******/////
    std::unique_ptr<TF1> fEtaPeak = std::make_unique<TF1>(Form("fEtaPeak%zu", iSlide), "gaus(0)+ pol2(3)", 0.40, 0.63); // gaus: [0]=A, [1]=mean, [2]=sigma 
    TH1D* hSignal = (TH1D*)hSame->Clone(Form("hSignal%zu", iSlide));
    hSignal->Add(hSignal, hMixScaled, 1., -1.);
    setHistoStandardSettings1D(hSignal);
    hSignal->Fit(fGaus.get(), "0M", "", 0.470, 0.570);        //TO DO: X WERTE VARIABLE DEFINIEREN!!
    fEtaPeak->SetParLimits(1, 0.95 * fGaus->GetParameter(1), 1.05 * fGaus->GetParameter(1));
    fEtaPeak->SetParLimits(2, 0.5 * fGaus->GetParameter(2), 1.05 * fGaus->GetParameter(2));
    hSignal->Fit(fEtaPeak.get(), "0M", "", 0.2, 0.8);
    hSignal->Write();
    fEtaPeak->Write();

    ////*****Calculating Yield*****/////
    ////***Data***/////
    double mean = fEtaPeak->GetParameter(1);
    double sigma = fEtaPeak->GetParameter(2);
    if(debug) out << "-------------CALULATING YIELD ----------- \n" << "mean of fEtaPeak: " << mean << "\n";
    if(debug) out << "sigma: " << sigma << "\n";

    // double xmin = mean - (3.0*sigma);
    // double xmax = mean + (3.0*sigma);

    double xmin = 0.53;
    double xmax = 0.57;

    int bmin = hSignal->GetXaxis()->FindBin(xmin);
    int bmax = hSignal->GetXaxis()->FindBin(xmax);

    double err = 0;

    double nEtas = hSignal->IntegralAndError(bmin, bmax, err);
    if(nEtas<0) {
      nEtas=0;
      err = 0;
    }
    if(debug) out << "nEtas: " << nEtas << "\n";
    hRawpT->SetBinContent(static_cast<Int_t>(iSlide) + 1, nEtas);
    hRawpT->SetBinError(static_cast<Int_t>(iSlide) + 1, err);

    if(isMC){
      ////***True***/////
      double errTrue = 0;
      double nTrueEtasPeak = hTrue->IntegralAndError(bmin, bmax, errTrue);
      // double nTrueEtas = hTrue->GetEntries();
      if(debug) out << "nTrueEtasPeak: " << nTrueEtasPeak << "\n";
      hTrueRawpT->SetBinContent(static_cast<Int_t>(iSlide) + 1, nTrueEtasPeak);
      hTrueRawpT->SetBinError(static_cast<Int_t>(iSlide) + 1, errTrue);
  
      ////***Generated Etas***/////
      double nGenEta = vecNGenEtaPerPt[iSlide];
      if(debug) out << "nGenEta: " << nGenEta << "\n";
      double nGenEtaError = vecNGenEtaPerPtError[iSlide];
      hGenRawpT->SetBinContent(static_cast<Int_t>(iSlide) + 1, nGenEta);
      hGenRawpT->SetBinError(static_cast<Int_t>(iSlide) + 1, nGenEtaError);
      
      double nGenEtaDalitz = nGenEta * 0.007;
      double errTrueBR = errTrue * 0.007;
      hGenRawpTDalitz->SetBinContent(static_cast<Int_t>(iSlide) + 1, nGenEtaDalitz);
      hGenRawpTDalitz->SetBinError(static_cast<Int_t>(iSlide) + 1, errTrueBR);
    }  
    
    ////*****Checks*****/////
    ////***Signal over Background***/////
    double errBack = 0 ;
    double background = hMixScaled->IntegralAndError(bmin, bmax, errBack);

    double sob = nEtas / background;
    double sobError = sqrt(((err/background) * (err/background)) + ((nEtas/(background*background)*errBack) * (nEtas/(background*background)*errBack)));    //delat f = sqrt( (delta x/ y)^2 + ( x/(y^2) * delta y )^2 )

    hSB->SetBinContent(static_cast<Int_t>(iSlide) + 1, sob);
    hSB->SetBinError(static_cast<Int_t>(iSlide) + 1, sobError);
    
    ////***Significance***/////
    if(isMC){
      double errSigni = 0;
      hSame->IntegralAndError(bmin, bmax, errSigni);

      if(errSigni > 0){
        double signi = nEtas/errSigni;
        hSigni->SetBinContent(static_cast<Int_t>(iSlide) + 1, signi);
      }
    }
    
  }
  hRawpT->Write();
  hSB->Write();
  if(isMC){
    hTrueRawpT->Write();
    hGenRawpT->Write();
    hGenRawpTDalitz->Write();
    hSigni->Write();
  }
}
