///////////////////////////////////////////////////////////////////////////////////
// Correcting Spectrum 
///////////////////////////////////////////////////////////////////////////////////
#include "helperPlotting.h"
#include <RtypesCore.h>

void correctingSpectra(const std::string& filenameRaws, const std::string& filenameEffi, const char* filenameCorrected = "corrFile.root", bool isMC = false){

    ////*****GETTING THE RAW YIELDS FROM EXTRACTION****/////
    TFile* fileFromRaws = safelyOpenRootfile(filenameRaws); 
    std::unique_ptr<TH1D> hRawpT(fileFromRaws->Get<TH1D>("hRawpT"));
    std::unique_ptr<TH1D> hNCollision(fileFromRaws->Get<TH1D>("hNCollision"));
    double_t nCollisions = hNCollision->GetBinContent(12);


    std::unique_ptr<TH1D> hTrueRawpT;
    std::unique_ptr<TH1D> hGenRawpT;
    std::unique_ptr<TH1D> hGenRawpTDalitz;

    if(isMC){
        auto* obj = fileFromRaws->Get<TH1D>("hTrueRawpT");
        hTrueRawpT.reset(static_cast<TH1D*>(obj->Clone("hTrueRawpTClone")));
        setHistoStandardSettings1D(hTrueRawpT.get());

        auto* obj2 = fileFromRaws->Get<TH1D>("hGenRawpT");
        hGenRawpT.reset(static_cast<TH1D*>(obj2->Clone("hGenRawpTClone")));
        setHistoStandardSettings1D(hGenRawpT.get());

        auto* obj3 = fileFromRaws->Get<TH1D>("hGenRawpTDalitz");
        hGenRawpTDalitz.reset(static_cast<TH1D*>(obj3->Clone("hGenRawpTDalitzCLone")));
        setHistoStandardSettings1D(hGenRawpTDalitz.get());
    }

    ////*****GETTING THE EFFIS****/////
    TFile* fileFromEffi = safelyOpenRootfile(filenameEffi); 
    std::unique_ptr<TH1D> hRecEffi(fileFromEffi->Get<TH1D>("hRecEffi"));
    std::unique_ptr<TH1D> hTrueEffi(fileFromEffi->Get<TH1D>("hTrueEffi"));
    std::unique_ptr<TH1D> hRecEffiDalitz(fileFromEffi->Get<TH1D>("hRecEffiDalitz"));
    std::unique_ptr<TH1D> hTrueEffiDalitz(fileFromEffi->Get<TH1D>("hTrueEffiDalitz"));

    ////*****CORRECTING RAW SPECTRA****/////
    std::unique_ptr<TFile> corrFile = std::unique_ptr<TFile>(TFile::Open(filenameCorrected, "RECREATE"));

    auto *hRecCorr = (TH1D*) hRawpT->Clone("hRecCorr");
    hRecCorr->Divide(hRecCorr, hRecEffi.get(), 1, 1, ""); // "B"??
    setHistoStandardSettings1D(hRecCorr);
    normalizeSepc(hRecCorr, nCollisions);
    hRecCorr->Write();

    if(isMC){
        auto *hTrueCorr = (TH1D*) hTrueRawpT->Clone("hTrueCorr");
        hTrueCorr->Divide(hTrueCorr, hTrueEffi.get(), 1, 1, ""); 
        setHistoStandardSettings1D(hTrueCorr);
        hTrueCorr->Write();
        
        auto *hTrueCorrNorm = (TH1D*) hTrueCorr->Clone("hTrueCorrNorm");
        normalizeSepc(hTrueCorrNorm, nCollisions);
        hTrueCorrNorm->Write();
    
        auto *hRecCorrDalitz = (TH1D*) hRawpT->Clone("hRecCorrDalitz");
        hRecCorrDalitz->Divide(hRecCorrDalitz, hRecEffiDalitz.get(), 1, 1, ""); 
        setHistoStandardSettings1D(hRecCorrDalitz);
        normalizeSepc(hRecCorrDalitz, nCollisions);
        hRecCorrDalitz->Write();
    
        auto *hTrueCorrDalitz = (TH1D*) hTrueRawpT->Clone("hTrueCorrDalitz");
        hTrueCorrDalitz->Divide(hTrueCorrDalitz, hTrueEffiDalitz.get(), 1, 1, ""); 
        setHistoStandardSettings1D(hTrueCorrDalitz);
        normalizeSepc(hTrueCorrDalitz, nCollisions);
        hTrueCorrDalitz->Write();
    }


    //TO DOO::
    // Normierungen noch einbauen!!
}
