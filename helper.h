#include <iostream>

// Header file to access Pythia 8 program elements.

/* Information about Plotting:
alle Messgroessen kursiv: p, N , m
nicht lursiv, _T, B!
*/

// ROOT, for histogramming und fitten
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TF1.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TGaxis.h" 
#include <TLine.h>
#include <TAttLine.h>
#include <TTree.h>
#include <TLatex.h>
#include <THnSparse.h>
#include "TColor.h"

//Canvas
#include <TCanvas.h>

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

//allgemein
#include <map>

//von mir: 
#include "TStyle.h"
#include "TMath.h"



//Standard settings
void SetCanvasStandardSettings(TCanvas *cCanv){
    cCanv->SetTopMargin(0.03);
    cCanv->SetBottomMargin(0.09);
    cCanv->SetRightMargin(0.025);
    cCanv->SetLeftMargin(0.09);
    cCanv->SetTickx();          //um "achsenstriche" oben und rechts an der Seite zu haben
    cCanv->SetTicky();
    cCanv->SetLogy(0);          // mit bool 
    cCanv->SetLogx(0);
    cCanv->SetLogz(0);
}

void SetSettingsCanvas_wRatio(TCanvas* c) {
    c->Divide(1,2);
    c->GetPad(1)->SetPad(0.,0.31,1,1);
    c->GetPad(2)->SetPad(0.,0.,1,0.31);     //x1,y1,x2,y2
  
    c->cd(1);
    gPad->SetTopMargin(0.07);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.0);
    gPad->SetTickx();
    gPad->SetTicky();
  
    c->cd(2);
    gPad->SetRightMargin(0.1);
    gPad->SetLeftMargin(0.1);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.4);
    gPad->SetTickx();
    gPad->SetTicky();

  }
  
void SetLegendSettings(TLegend* leg){
    leg->SetTextSize(27);
    leg->SetTextFont(43);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetLineColor(0);
    leg->SetMargin(0.15);
}
void SetGlobalStandard(void)
{
    gStyle->SetOptStat(0);           //~the box is not shown
    gStyle->SetCanvasColor(0);       //~the canvas is white (=0)
    gStyle->SetPadColor(0);          //~sets the pads, a part of the canvas, to the color white
    gStyle->SetCanvasBorderMode(0);  //~sets the border of the canvas to white 
    gStyle->SetPadBorderMode(0);     //~the thickness of the border is 0
    gStyle->SetNumberContours(256);
}

void SetHeaderSettings(TLatex* head){
    head->SetTextSize(27);
    head->SetTextFont(43);
}

void SetHistoStandardSettings1D(TH1* histo, double XOffset = 1.2, double YOffset = 1.1){
    //Title
    histo->SetTitle("");

    histo->GetXaxis()->SetTitleOffset(XOffset);
    histo->GetYaxis()->SetTitleOffset(YOffset);
    // histo->GetXaxis()->SetTitleSize(35);      //achsen beschriftung
    // histo->GetYaxis()->SetTitleSize(35);
    // histo->GetYaxis()->SetTitleFont(43);
    // histo->GetXaxis()->SetTitleFont(43);
    histo->GetXaxis()->SetTitleSize(25); 
    histo->GetYaxis()->SetTitleSize(25);  
    histo->GetXaxis()->SetTitleFont(63); 
    histo->GetYaxis()->SetTitleFont(63); 

    //Label
    histo->GetXaxis()->SetLabelSize(23);      //label = zahlen an der achse
    histo->GetYaxis()->SetLabelSize(23);
    histo->GetXaxis()->SetLabelFont(43);      //evtl 42, aber 43 ist okay
    histo->GetYaxis()->SetLabelFont(43);

    //Marker
    histo->SetMarkerStyle(20);                //zeichen, das den Punkt anzeigt. 20 ist kreis, gibt auch rechteck usw vgl. internet
    histo->SetMarkerSize(1);                //By default auf 1, kann beliebig verändert werden, also auch auf 1.25 usw
    histo->SetLineWidth(3);                   
    histo->SetLineColor(kBlack);
    histo->SetMarkerColor(kBlack);
}

void SetHistoStandardSettings2D(TH2* histo, double XOffset = 1.2, double YOffset = 1.5, double ZOffset = 1.3){
    //Title
    histo->GetXaxis()->SetTitleOffset(XOffset);
    histo->GetYaxis()->SetTitleOffset(YOffset);
    histo->GetZaxis()->SetTitleOffset(ZOffset);
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

    //Label
    histo->GetXaxis()->SetLabelSize(23);    //originally 45
    histo->GetYaxis()->SetLabelSize(23);    //originally 45
    histo->GetZaxis()->SetLabelSize(23); 
    histo->GetXaxis()->SetLabelFont(43);
    histo->GetYaxis()->SetLabelFont(43);
    histo->GetZaxis()->SetLabelFont(43);

}


//Functions
// Function to open a ROOT file without changing the gDirectory
TFile* SafelyOpenRootfile(const std::string filename){
    // if (!filename || filename->IsZombie()) {
    //     std::cerr << "MC file konnte nicht geöffnet werden!\n";
    // }
  /// Opens a rootfile without affecting the active path, which otherwise would point into the file, often causing trouble.
  //save current path before opening rootfile.
  TString sPath = gDirectory->GetPath();

  TFile* ffile = 0x0;
  // check if file is already open.
  if ( gROOT->GetListOfFiles()->FindObject(filename.data()) ) {
    ffile = gROOT->GetFile(filename.data()); // avoid to open same file twice
  }
  if (ffile && ffile->IsOpen()) return ffile;

  ffile = TFile::Open(filename.data()); // gives root error and returns 0x0 on fail.
  // if (!ffile) printf(Form("SafelyOpenRootfile(): file '%s' not found.", filename.data()));

  // change to previous path again, so that it will be possible to close the file later without crash.
  // otherwise heap based objects will be created in memory that will be freed when the file is closed.
  gDirectory->Cd(sPath.Data());

  // alternatively one can do hist->SetDirectory(0); // according to Dario (Analysis Tutorial 26.06.2015)
  // but this seems not to work if the class object (which owns the hist) was created in the file path.

  return ffile;

  
}

//sortingVectors
void SortingVectorN3(int NSparse, TH1D* histo, std::vector<TH1D*>& vec1, std::vector<TH1D*>& vec2, std::vector<TH1D*>& vec3){
    if( NSparse == 0){
        vec1.push_back(histo);
    }
    if( NSparse == 1){
        vec2.push_back(histo);
    }
    if( NSparse == 2){
        vec3.push_back(histo);
    }
    if(NSparse > 2) {
        std::cout << " NSparse zu groß!" << std::endl;
    }
}

//Making Ratio
TH1D* MakingRatio( TH1* histo, TH1* histo_divide, const char* name){
    TH1D *hist_ratio = (TH1D*)histo->Clone(name);
    hist_ratio->Divide(hist_ratio, histo_divide, 1. , 1., "B");   // B: berücksigtigt korrelationen
  
    // TH1D *hist_ratio = (TH1D*)histo_clo->Clone(name);
    SetHistoStandardSettings1D(hist_ratio, 1.3, 1.4);
    return hist_ratio;
}

//Fitting
//Crystal ball function for signal +linear background, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBallBck(Double_t *x,Double_t *par) { //Funktioniert mit den Paramtern unter okay gut

    Double_t t = (x[0]-par[1])/par[2];
    if (par[4] < 0) t = -t;

    Double_t absAlpha = fabs((Double_t)par[4]);

    if (t >= -absAlpha) {
        return par[0]*exp(-0.5*t*t)+par[5]+par[6]*x[0];
    } else {
        Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
        Double_t b= par[3]/absAlpha - absAlpha;

        return par[0]*(a/TMath::Power(b - t, par[3]))+par[5]+par[6]*x[0];
    }
}
//for cleaning the code, regards to CrystalBallBck
void SetFitParas(TF1* fit){
    // double mesonAmplitude = h_signal->GetMaximum(); //verschiebt maximum
    double expectmass = 0.126;
    double expectedwidth = 0.007;

    fit->SetParLimits(1, 0.5, 6.0);
    fit->SetParLimits(2, 1.1, 10.0);
    fit->SetParLimits(3, 0.003, 0.020);
    fit->SetParLimits(4, 0.125, 0.145);
    fit->SetParLimits(6, -100.0, -1e-3);
    fit->SetParameter(3,4.5);   //original: 2.
    fit->SetParameter(4, 0.7);  //original: 0.7  //wenn groß, dann quasi kein tail nach links, wenn klein, dann großen tail nach links
    fit->SetParameter(5,15);
    fit->SetParameter(6,1.);
    fit->SetLineWidth(3); 
    fit->SetNpx(10000);
    fit->SetLineColor(kBlue);
}


//second function (from Steffi)
static Double_t CrystalBallLow(Double_t x, Double_t alpha, Double_t n,
    Double_t sigma, Double_t mean) {
const Double_t t = (x - mean) / sigma;
const Double_t absA = std::fabs(alpha);
if (t > -absA)
return std::exp(-0.5 * t * t);
const Double_t a = std::pow(n / absA, n) * std::exp(-0.5 * absA * absA);
const Double_t b = n / absA - absA;
return a / std::pow(b - t, n);
}

static Double_t CBExpoModel(Double_t *xx, Double_t *par) {
    const Double_t x = xx[0];
    return par[0] * CrystalBallLow(x, par[1], par[2], par[3], par[4]) + std::exp(par[5] + par[6] * x);
}

static TF1 *FitPi0Shape_Save(TH1 *hSig, double fitMin, double fitMax) {     //funktioniert gut
    if (!hSig || hSig->GetNbinsX() == 0) {
        return nullptr;
    }

    TF1 fitFunc(Form("fit_%s_tmp", hSig->GetName()), CBExpoModel, fitMin, fitMax, 7);

    const double ymax = std::max(1.0, hSig->GetMaximum());
    fitFunc.SetParameters(
        ymax,                                   // par[0] Amplitude CB
        1.5,                                    // par[1] alpha
        3.0,                                    // par[2] n
        0.008,                                  // par[3] sigma
        0.135,                                  // par[4] mean (pi0)
        std::log(std::max(1.0, ymax / 2.0)),    // par[5] log-Amplitude Untergrund
        -8.0                                    // par[6] slope exp.
    );

    fitFunc.SetParLimits(1, 0.5, 6.0);
    fitFunc.SetParLimits(2, 1.1, 10.0);
    fitFunc.SetParLimits(3, 0.003, 0.020);
    fitFunc.SetParLimits(4, 0.125, 0.145);
    fitFunc.SetParLimits(6, -100.0, -1e-3);
    fitFunc.SetNpx(10000);

    hSig->Fit(&fitFunc, "RQ0M");

    TF1 *fcl = (TF1 *)fitFunc.Clone(Form("fit_%s", hSig->GetName()));
    if (fcl) {
        fcl->SetRange(fitMin, fitMax);
        hSig->GetListOfFunctions()->Add(fcl);
    }

    return fcl;
}

void DrawingHeaderStandardLines(double_t x, double_t y1, double_t y2, double_t y3){
    std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>(); 		
	SetHeaderSettings(lHeader.get());

    lHeader->DrawLatexNDC(x,y1, "ALICE work in progress");
    lHeader->DrawLatexNDC(x, y2, "LHC23zo,zp");
    lHeader->DrawLatexNDC(x, y3, "pp, #sqrt{#it{s}} = 13.6 TeV, B = 0.2 T");
}

void DrawingHeaderStandardLinesMC(double_t x, double_t y1, double_t y2, double_t y3){
    std::unique_ptr<TLatex> lHeader = std::make_unique<TLatex>(); 		
	SetHeaderSettings(lHeader.get());

    lHeader->DrawLatexNDC(x,y1, "ALICE work in progress");
    lHeader->DrawLatexNDC(x, y2, "MC, LHC25d2");
    lHeader->DrawLatexNDC(x, y3, "pp, #sqrt{#it{s}} = 13.6 TeV, B = 0.2 T");
}

void SavingPdfPng(const char *filename="plot") {
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas", "", 800, 800);
    SetCanvasStandardSettings(canvas.get());
    canvas->SaveAs(Form("Plots_zozp_fsp/%s.png", filename));
    canvas->SaveAs(Form("/Users/lauragans-bartl/Downloads/%s.pdf", filename));
}