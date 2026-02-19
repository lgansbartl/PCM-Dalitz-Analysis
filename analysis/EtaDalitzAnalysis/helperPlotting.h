// Header file to access Pythia 8 program elements.

/* Information about Plotting:
alle Messgroessen kursiv: p, N , m
nicht lursiv, _T, B!
*/
#ifndef HELPER_PLOTTING_H
#define HELPER_PLOTTING_H

// ROOT, for histogramming und fitten
#include "TColor.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"
#include <TAttLine.h>
#include <THnSparse.h>
#include <TLatex.h>
#include <TLine.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include "stdarg.h"
#include "TMath.h"

// Canvas
#include <TCanvas.h>

// ROOT, for interactive graphics.
#include "TApplication.h"
#include "TVirtualPad.h"

// ROOT, for saving file.
#include "TFile.h"

// allgemein
#include <iostream>
#include <map>

// von mir:
#include "TMath.h"
#include "TStyle.h"

// Standard settings
void setCanvasStandardSettings(TCanvas* cCanv)
{
  cCanv->SetTopMargin(static_cast<float>(0.05));
  cCanv->SetBottomMargin(static_cast<float>(0.09));
  cCanv->SetRightMargin(static_cast<float>(0.025));
  cCanv->SetLeftMargin(static_cast<float>(0.09));
  cCanv->SetTickx(); // um "achsenstriche" oben und rechts an der Seite zu haben
  cCanv->SetTicky();
  cCanv->SetLogy(0); // mit bool
  cCanv->SetLogx(0);
  cCanv->SetLogz(0);
}

void setSettingsCanvaswRatio(TCanvas* canvas)
{
  canvas->Divide(1, 2);
  canvas->GetPad(1)->SetPad(0., 0.31, 1, 1);
  canvas->GetPad(2)->SetPad(0., 0., 1, 0.31); // x1,y1,x2,y2

  canvas->cd(1);
  gPad->SetTopMargin(static_cast<float>(0.07));
  gPad->SetRightMargin(static_cast<float>(0.1));
  gPad->SetLeftMargin(static_cast<float>(0.1));
  gPad->SetBottomMargin(0.0);
  gPad->SetTickx();
  gPad->SetTicky();

  canvas->cd(2);
  gPad->SetRightMargin(static_cast<float>(0.1));
  gPad->SetLeftMargin(static_cast<float>(0.1));
  gPad->SetTopMargin(static_cast<float>(0.0));
  gPad->SetBottomMargin(static_cast<float>(0.4));
  gPad->SetTickx();
  gPad->SetTicky();
}

void setLegendSettings(TLegend* leg)
{
  leg->SetTextSize(27);
  leg->SetTextFont(43);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  leg->SetMargin(static_cast<float>(0.15));
}
void setGlobalStandard()
{
  gStyle->SetOptStat(0);          //~the box is not shown
  gStyle->SetCanvasColor(0);      //~the canvas is white (=0)
  gStyle->SetPadColor(0);         //~sets the pads, a part of the canvas, to the color white
  gStyle->SetCanvasBorderMode(0); //~sets the border of the canvas to white
  gStyle->SetPadBorderMode(0);    //~the thickness of the border is 0
  gStyle->SetNumberContours(256);
}

void setHeaderSettings(TLatex* head)
{
  head->SetTextSize(27);
  head->SetTextFont(43);
}

void setHistoStandardSettings1D(TH1* histo, double xOffset = 1.2, double yOffset = 1.1)
{
  // Title
  histo->SetTitle("");

  histo->GetXaxis()->SetTitleOffset(static_cast<float>(xOffset));
  histo->GetYaxis()->SetTitleOffset(static_cast<float>(yOffset));
  // histo->GetXaxis()->SetTitleSize(35);      //achsen beschriftung
  // histo->GetYaxis()->SetTitleSize(35);
  // histo->GetYaxis()->SetTitleFont(43);
  // histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetTitleFont(63);
  histo->GetYaxis()->SetTitleFont(63);

  // Label
  histo->GetXaxis()->SetLabelSize(23); // label = zahlen an der achse
  histo->GetYaxis()->SetLabelSize(23);
  histo->GetXaxis()->SetLabelFont(43); // evtl 42, aber 43 ist okay
  histo->GetYaxis()->SetLabelFont(43);

  // Marker
  histo->SetMarkerStyle(20); // zeichen, das den Punkt anzeigt. 20 ist kreis, gibt auch rechteck usw vgl. internet
  histo->SetMarkerSize(1);   // By default auf 1, kann beliebig verändert werden, also auch auf 1.25 usw
  histo->SetLineWidth(3);
  histo->SetLineColor(kBlack);
  histo->SetMarkerColor(kBlack);
}

void setHistoStandardSettings2D(TH2* histo, double xOffset = 1.2, double yOffset = 1.5, double zOffset = 1.3)
{
  // Title
  histo->GetXaxis()->SetTitleOffset(static_cast<float>(xOffset));
  histo->GetYaxis()->SetTitleOffset(static_cast<float>(yOffset));
  histo->GetZaxis()->SetTitleOffset(static_cast<float>(zOffset));
  // histo->GetXaxis()->SetTitleSize(35);    //originally 45
  // histo->GetYaxis()->SetTitleSize(35);    //originally 45
  // histo->GetZaxis()->SetTitleSize(35);
  // histo->GetXaxis()->SetTitleFont(43);
  // histo->GetYaxis()->SetTitleFont(43);
  // histo->GetZaxis()->SetTitleFont(43);

  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetZaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetTitleFont(63);
  histo->GetYaxis()->SetTitleFont(63);
  histo->GetZaxis()->SetTitleFont(63);

  // Label
  histo->GetXaxis()->SetLabelSize(23); // originally 45
  histo->GetYaxis()->SetLabelSize(23); // originally 45
  histo->GetZaxis()->SetLabelSize(23);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetZaxis()->SetLabelFont(43);
}

// Functions
//  Function to open a ROOT file without changing the gDirectory
TFile* safelyOpenRootfile(const std::string &filename)
{
  // if (!filename || filename->IsZombie()) {
  //     std::cerr << "MC file konnte nicht geöffnet werden!\n";
  // }
  /// Opens a rootfile without affecting the active path, which otherwise would point into the file, often causing trouble.
  // save current path before opening rootfile.
  TString sPath = gDirectory->GetPath();

  TFile* ffile = nullptr;
  // check if file is already open.
  if (gROOT->GetListOfFiles()->FindObject(filename.data())) {
    ffile = gROOT->GetFile(filename.data()); // avoid to open same file twice
  }
  if (ffile != nullptr && ffile->IsOpen()) {
    ffile->cd();
    return ffile;
  }
 
  ffile = TFile::Open(filename.data()); // gives root error and returns 0x0 on fail.
  // if (!ffile) printf(Form("SafelyOpenRootfile(): file '%s' not found.", filename.data()));

  // change to previous path again, so that it will be possible to close the file later without crash.
  // otherwise heap based objects will be created in memory that will be freed when the file is closed.
  gDirectory->Cd(sPath.Data());

  // alternatively one can do hist->SetDirectory(0); // according to Dario (Analysis Tutorial 26.06.2015)
  // but this seems not to work if the class object (which owns the hist) was created in the file path.

  return ffile;
}

// sortingVectors
void sortingVectorN3(int nSparse, TH1D* histo, std::vector<TH1D*>& vec1, std::vector<TH1D*>& vec2, std::vector<TH1D*>& vec3)
{
  if (nSparse == 0) {
    vec1.push_back(histo);
  }
  if (nSparse == 1) {
    vec2.push_back(histo);
  }
  if (nSparse == 2) {
    vec3.push_back(histo);
  }
  if (nSparse > 2) {
    std::cout << " nSparse zu groß \n!";
  }
}

// Making Ratio ->evtl in struct umbauen?
TH1D* makingRatio(TH1* histo, TH1* histo_divide, const char* name)
{
  TH1D* hist_ratio = reinterpret_cast<TH1D*>(histo->Clone(name));
  hist_ratio->Divide(hist_ratio, histo_divide, 1., 1., "B"); // B: berücksigtigt korrelationen

  // TH1D *hist_ratio = (TH1D*)histo_clo->Clone(name);
  setHistoStandardSettings1D(hist_ratio, 1.3, 1.4);
  return hist_ratio;
}

// Fitting
// Crystal ball function for signal +linear background, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t crystalBallBck(Double_t* xIn, Double_t* par)
{ // Funktioniert mit den Paramtern unter okay gut

  Double_t var = (xIn[0] - par[1]) / par[2];

  if (par[4] < 0){
    var = -var;
  }

  Double_t absAlpha = fabs(par[4]);

  if (var >= -absAlpha) {
    return (par[0] * exp(-0.5 * var * var)) + par[5] + (par[6] * xIn[0]);
  } 
  Double_t aVar = TMath::Power(par[3] / absAlpha, par[3]) * exp(-0.5 * absAlpha * absAlpha);
  Double_t bVar = (par[3] / absAlpha) - absAlpha;

  return (par[0] * (aVar / TMath::Power(bVar - var, par[3]))) + par[5] + (par[6] * xIn[0]);
}
// for cleaning the code, regards to CrystalBallBck
void setFitParas(TF1* fit)
{
  // double mesonAmplitude = h_signal->GetMaximum(); //verschiebt maximum
  // double expectmass = 0.126;
  // double expectedwidth = 0.007;

  fit->SetParLimits(1, 0.5, 6.0);
  fit->SetParLimits(2, 1.1, 10.0);
  fit->SetParLimits(3, 0.003, 0.020);
  fit->SetParLimits(4, 0.125, 0.145);
  fit->SetParLimits(6, -100.0, -1e-3);
  fit->SetParameter(3, 4.5); // original: 2.
  fit->SetParameter(4, 0.7); // original: 0.7  //wenn groß, dann quasi kein tail nach links, wenn klein, dann großen tail nach links
  fit->SetParameter(5, 15);
  fit->SetParameter(6, 1.);
  fit->SetLineWidth(3);
  fit->SetNpx(10000);
  fit->SetLineColor(kBlue);
}

// second function (from Steffi)
static Double_t crystalBallLow(Double_t xVar, Double_t alpha, Double_t nVar, Double_t sigma, Double_t mean)
{
  const Double_t tVar = (xVar - mean) / sigma;
  const Double_t absA = std::fabs(alpha);
  if (tVar > -absA) {
    return std::exp(-0.5 * tVar * tVar);
  }
  const Double_t aVar = std::pow(nVar / absA, nVar) * std::exp(-0.5 * absA * absA);
  const Double_t bVar = (nVar / absA) - absA;
  return aVar / std::pow(bVar - tVar, nVar);
}

static Double_t cbExpoModel(const Double_t* var, const Double_t* par)
{
  const Double_t var0 = var[0];
  return (par[0] * crystalBallLow(var0, par[1], par[2], par[3], par[4])) + std::exp(par[5] + (par[6] * var0));
}

static TF1* fitPi0ShapeSave(TH1* hSig, double fitMin, double fitMax)
{ // funktioniert gut
  if (hSig == nullptr|| hSig->GetNbinsX() == 0) {
    return nullptr;
  }

  TF1 fitFunc(Form("fit_%s_tmp", hSig->GetName()), cbExpoModel, fitMin, fitMax, 7);

  const double ymax = std::max(1.0, hSig->GetMaximum());
  fitFunc.SetParameters(
    ymax,                                // par[0] Amplitude CB
    1.5,                                 // par[1] alpha
    3.0,                                 // par[2] n
    0.008,                               // par[3] sigma
    0.135,                               // par[4] mean (pi0)
    std::log(std::max(1.0, ymax / 2.0)), // par[5] log-Amplitude Untergrund
    -8.0                                 // par[6] slope exp.
  );

  fitFunc.SetParLimits(1, 0.5, 6.0);
  fitFunc.SetParLimits(2, 1.1, 10.0);
  fitFunc.SetParLimits(3, 0.003, 0.020);
  fitFunc.SetParLimits(4, 0.125, 0.145);
  fitFunc.SetParLimits(6, -100.0, -1e-3);
  fitFunc.SetNpx(10000);

  hSig->Fit(&fitFunc, "RQ0M");

  TF1* fcl = (TF1*)fitFunc.Clone(Form("fit_%s", hSig->GetName()));
  if (fcl != nullptr) {
    fcl->SetRange(fitMin, fitMax);
    hSig->GetListOfFunctions()->Add(fcl);
  }

  return fcl;
}

void drawingHeaderStandardLines(double_t xVar, double_t y1Var, const char* dataset, const char* cllisionAndbField, const char* particleDecay)
{
  std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>();
  setHeaderSettings(lHeader.get());

  double_t dy = 0.05;
  double_t y2Var = y1Var - dy;
  double_t y3Var = y1Var - (2*dy);
  double_t y4Var = y1Var - (3*dy);

  lHeader->DrawLatexNDC(xVar, y1Var, "ALICE work in progress");
  lHeader->DrawLatexNDC(xVar, y2Var, dataset);
  lHeader->DrawLatexNDC(xVar, y3Var, cllisionAndbField);
  lHeader->DrawLatexNDC(xVar, y4Var, particleDecay);

}

void drawingMarker(TH1* histo, Color_t color = kBlack, Style_t style = 20, Size_t size = 1.){
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(size);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

void savingPdfPng(const char* filename = "plot")
{
  std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
  setCanvasStandardSettings(canvas.get());
  canvas->SaveAs(Form("Plots_zozp_fsp/%s.png", filename));
  canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/%s.pdf", filename));
}

void drawLine(double x1l, double y1l, double x2l, double y2l, short width, short color){
  std::unique_ptr<TLine> line = std::make_unique<TLine> (x1l , y1l, x2l, y2l);
  line->SetLineWidth(width);
  line->SetLineColor(color);
  line->Draw("same");
}

void normalizeSepc(TH1* histo, double_t nCollisions = 1){
  histo->Scale(1, "width");   //Binbreite
  histo->Scale(1/(2* TMath::Pi()));      //Phi
  histo->Scale(1/nCollisions);            //Kollisionen
  histo->Scale(1/0.007);              //BR            
  histo->Scale(1/1.6);              //eta: |eta| < 0.8
  histo->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{1}{2#pi #it{p}_{T}} #frac{d^{2} #it{N}}{d#it{p}_{T}dy} #frac{1}{BR_{Dalitz}}((GeV/#it{c})^{-2})");

  //Scaling pT
  std::vector<Int_t> vecBins;
  for(int iBin = 1; iBin <= histo->GetNbinsX(); iBin++){
      double binCenter = histo->GetBinCenter(iBin); 
      if(binCenter == 0) continue;

      double xLow = histo->GetXaxis()->GetBinLowEdge(iBin);
      vecBins.push_back(static_cast<Int_t>(xLow));

      double binContent = histo->GetBinContent(iBin);
      double binError = histo->GetBinError(iBin);

      histo->SetBinContent(iBin, binContent / binCenter);
      histo->SetBinError(iBin, binError / binCenter);
  }
}


#endif // HELPER_PLOTTING_H
