
//******************INCLUDINGS******************//
#include "ExtractSignal.hxx"
#include "CalculateEffi.hxx"
#include "CorrectSpectrum.h"

int main(){
  
    //******************PATHS TO FILES******************//
    //***low  field 2023***//
    // const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC23zozp.root";        //merged

    //***Nominal fild 2024***//
    const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24.root";
    const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24f4dR.root";
    const char* filenameMCTrue = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24f4dT.root";


    //******************SIGNAL EXTRACTION******************//
    extractSignal(filename, "", "extrFilesData.root");

    bool isMC =true;
    extractSignal(filenameMCrec, filenameMCTrue, "extrFilesMC.root", isMC);

    // //******************CALCULATING EFFI******************//
    const char* extrFilesMC = "extrFilesMC.root";
    calculateEffi(extrFilesMC, "effiFiles.root");

    // //******************CORRECTIN SPECTRA******************//
    //files with results from previous steps
    const char* rawYieldFilesData = "extrFilesData.root";
    const char* rawYieldFilesMC = "extrFilesMC.root";
    const char* effiFiles = "effiFiles.root";  
    
    //names for files from this step
    const char* corrFilesData = "corrFilesData.root";
    const char* corrFilesMC = "corrFilesMC.root";

    correctingSpectra(rawYieldFilesData, effiFiles, corrFilesData);
    correctingSpectra(rawYieldFilesMC, effiFiles, corrFilesMC, isMC);
    
}
