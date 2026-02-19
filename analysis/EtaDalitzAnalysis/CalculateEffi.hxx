///////////////////////////////////////////////////////////////////////////////////
// Calculating Raw Yield and Efficiency 
///////////////////////////////////////////////////////////////////////////////////
#include "helperPlotting.h"

void calculateEffi(const std::string& filename, const char* filenameEffi){
    ////**********GETTING FILES**********/////
    TFile* fileFromExtraction = safelyOpenRootfile(filename); 
    std::unique_ptr<TH1D> hRawpT(fileFromExtraction->Get<TH1D>("hRawpT"));
    std::unique_ptr<TH1D> hTrueRawpT(fileFromExtraction->Get<TH1D>("hTrueRawpT"));
    std::unique_ptr<TH1D> hGenRawpT(fileFromExtraction->Get<TH1D>("hGenRawpT"));
    std::unique_ptr<TH1D> hGenRawpTDalitz(fileFromExtraction->Get<TH1D>("hGenRawpTDalitz"));

    ////**********GDEFINITIONS**********/////
    std::unique_ptr<TFile> effiFile = std::unique_ptr<TFile>(TFile::Open(filenameEffi, "RECREATE"));
    hRawpT->Write();
    hTrueRawpT->Write();
    hGenRawpT->Write();
    hGenRawpTDalitz->Write();

    
   
    ////**********Calculating EffICIENCY **********/////
    // Commmand to error propagation: the two hispts ar ekorrelated but not fully. When using the Divide() function, root expects the histos to be not correalted. 
    // Since they are, the errors are overestimated. The option B propagates the error binominally, wich is correct for fully correalted histos.
    //Since they are not fully correlated, using this option leads to underestimated errors
    // -> both options are not ideal!!!

    auto *hRecEffi = (TH1D*) hRawpT->Clone("hRecEffi");
    hRecEffi->Divide(hRecEffi, hGenRawpT.get(), 1, 1, ""); // "B"??
    setHistoStandardSettings1D(hRecEffi);
    hRecEffi->Write();

    auto *hTrueEffi = (TH1D*) hTrueRawpT->Clone("hTrueEffi");
    hTrueEffi->Divide(hTrueEffi, hGenRawpT.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hTrueEffi);
    hTrueEffi->Write();

    auto *hRecEffiDalitz = (TH1D*) hRawpT->Clone("hRecEffiDalitz");
    hRecEffiDalitz->Divide(hRecEffiDalitz, hGenRawpTDalitz.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hRecEffiDalitz);
    hRecEffiDalitz->Write();

    auto *hTrueEffiDalitz = (TH1D*) hTrueRawpT->Clone("hTrueEffiDalitz");
    hTrueEffiDalitz->Divide(hTrueEffiDalitz, hGenRawpTDalitz.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hTrueEffiDalitz);
    hTrueEffiDalitz->Write();

}

