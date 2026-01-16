void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary, TString functionname)
{

  fFileErrLog << "Start Fitting spectra with CB fit" << endl;
  fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0], fMesonMassRange[1]);
  Double_t mesonAmplitude = fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
  Double_t mesonAmplitudeMin = mesonAmplitude * 50. / 100.;
  Double_t mesonAmplitudeMax = mesonAmplitude * 115. / 100.;

  fFitReco = new TF1(functionname, CrystalBallBck, fMesonFitRange[0], fMesonFitRange[1], 7);

  fFitGausExp = new TF1("crystal", CrystalBall, fMesonFitRange[0], fMesonFitRange[1], 5);

  fFitLinearBck = new TF1("Linear", "[0]+[1]*x", fMesonFitRange[0], fMesonFitRange[1]);

  fFitReco->SetParameter(0, mesonAmplitude);
  fFitReco->SetParameter(1, fMesonMassExpect);
  fFitReco->SetParameter(2, fMesonWidthExpect);
  fFitReco->SetParameter(3, 2.);
  fFitReco->SetParameter(4, 0.7);
  fFitReco->SetParameter(5, 0.);
  fFitReco->SetParameter(6, 1.);

  if (!vary) {
    fFitReco->FixParameter(3, fCBn);
    fFitReco->FixParameter(4, fCBAlpha);
  }

  fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);
  fFitReco->SetParLimits(1, fMesonMassRange[0], fMesonMassRange[1]);
  fFitReco->SetParLimits(2, fMesonWidthRange[0], fMesonWidthRange[1]);

  fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco, "QR0");
  fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco, "QRE0");

  fFitReco->SetLineColor(3);
  fFitReco->SetLineWidth(1);
  fFitReco->SetLineStyle(1);

  if (vary) {
    fCBAlpha = fFitReco->GetParameter(4);
    fCBn = fFitReco->GetParError(3);
  }

  fFitGausExp->SetParameter(0, fFitReco->GetParameter(0));
  fFitGausExp->SetParameter(1, fFitReco->GetParameter(1));
  fFitGausExp->SetParameter(2, fFitReco->GetParameter(2));
  fFitGausExp->SetParameter(3, fFitReco->GetParameter(3));
  fFitGausExp->SetParameter(4, fFitReco->GetParameter(4));

  fFitGausExp->SetParError(0, fFitReco->GetParError(0));
  fFitGausExp->SetParError(1, fFitReco->GetParError(1));
  fFitGausExp->SetParError(2, fFitReco->GetParError(2));
  fFitGausExp->SetParError(3, fFitReco->GetParError(3));
  fFitGausExp->SetParError(4, fFitReco->GetParError(4));

  fFitLinearBck->SetParameter(0, fFitReco->GetParameter(5));
  fFitLinearBck->SetParameter(1, fFitReco->GetParameter(6));

  fFitLinearBck->SetParError(0, fFitReco->GetParError(5));
  fFitLinearBck->SetParError(1, fFitReco->GetParError(6));

  Int_t binCenterStart;
  Double_t startBinEdge;
  Int_t binCenterEnd;
  Double_t endBinEdge;

  TVirtualFitter* fitter = TVirtualFitter::GetFitter();

  if (TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0) {
    binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[0] - (fMesonMassExpect - fFitReco->GetParameter(1)));
    startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart) - 0.5 * fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntDeltaRangeFit[1] - (fMesonMassExpect - fFitReco->GetParameter(1)));
    endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd) + 0.5 * fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

    Int_t nFreePar = fFitReco->GetNumberFreeParameters();
    double* covMatrix = fitter->GetCovarianceMatrix();

    Float_t intLinearBack = fFitLinearBck->GetParameter(0) * (endBinEdge - startBinEdge) +
                            0.5 * fFitLinearBck->GetParameter(1) * (endBinEdge * endBinEdge - startBinEdge * startBinEdge);

    Float_t errorLinearBck = pow((pow((endBinEdge - startBinEdge) * fFitReco->GetParError(5), 2) + pow(0.5 * (endBinEdge * endBinEdge - startBinEdge * startBinEdge) * fFitReco->GetParError(6), 2) + 2 * covMatrix[nFreePar * nFreePar - 2] * (endBinEdge - startBinEdge) * 0.5 * (endBinEdge * endBinEdge - startBinEdge * startBinEdge)), 0.5);

    fFileDataLog << "Parameter for bin " << ptBin << endl;
    fFileDataLog << "CrystalBall: \t" << fFitReco->GetParameter(0) << "+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1) << "+-" << fFitReco->GetParError(1) << "\t " << fFitReco->GetParameter(2) << "+-" << fFitReco->GetParError(2) << "\t " << fFitReco->GetParameter(3) << "+-" << fFitReco->GetParError(3) << "\t " << fFitReco->GetParameter(4) << "+-" << fFitReco->GetParError(4) << endl;
    fFileDataLog << "Linear: \t" << fFitReco->GetParameter(5) << "+-" << fFitReco->GetParError(5) << "\t " << fFitReco->GetParameter(6) << "+-" << fFitReco->GetParError(6) << endl;

    fIntLinearBck = intLinearBack / fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    fIntLinearBckError = errorLinearBck / fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
  } else {
    fFileErrLog << "Fitting failed in " << ptBin << " with status::" << gMinuit->fCstatu.Data() << "why failed?" << endl
                << endl;
  }
  fFitReco->DrawCopy("same");
}
