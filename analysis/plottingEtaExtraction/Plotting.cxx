#include "PlotManager.h"
#include "Logging.h"
#include <string>

using std::string;
using namespace SciRooPlot;
void DefineInputIdentifiers(PlotManager& plotManager);
void DefinePlotTemplates(PlotManager& plotManager);
void DefinePlots(PlotManager& plotManager);


//****************************************************************************************
/*
Global Definitions
*/
string dataset = "LHC24";

//****************************************************************************************

int main(int argc, char* argv[])
{
  PlotManager plotManager;
  DefineInputIdentifiers(plotManager);
  DefinePlotTemplates(plotManager);
  DefinePlots(plotManager);
  plotManager.SaveProject("plottingEtaExtraction");
  return 0;
}

//****************************************************************************************
/**
 * Defines the input identifiers.
 */
//****************************************************************************************
void DefineInputIdentifiers(PlotManager& plotManager)
{
  const string inputFolder = "/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/build/";
  plotManager.AddInputDataFiles("extrFileData", {inputFolder + "extrFilesData.root"});
  plotManager.AddInputDataFiles("extrFileMC", {inputFolder + "extrFilesMC.root"});
  plotManager.AddInputDataFiles("effiFiles", {inputFolder + "effiFiles.root"});
  plotManager.AddInputDataFiles("corrFilesData", {inputFolder + "corrFilesData.root"});
  plotManager.AddInputDataFiles("corrFilesMC", {inputFolder + "corrFilesMC.root"});

}

//****************************************************************************************
/**
 * Defines the plot templates.
 */
//****************************************************************************************
void DefinePlotTemplates(PlotManager& plotManager)
{
  {
    auto plotTemplate = PlotManager::GetPlotTemplate("1d");
    plotTemplate[0].SetMargins(0.07, 0.14, 0.12, 0.07);
    plotManager.AddPlotTemplate(plotTemplate);
  }
  {
    auto plotTemplate = PlotManager::GetPlotTemplate("1d_ratio");
    plotManager.AddPlotTemplate(plotTemplate);
  }
  {
    auto plotTemplate = PlotManager::GetPlotTemplate("2d");
    plotTemplate[0].SetMargins(0.07, 0.14, 0.12, 0.07);
    plotManager.AddPlotTemplate(plotTemplate);
  }
  {  // -----------------------------------------------------------------------
    Plot plotTemplate("lMyTemplate", "PLOT_TEMPLATES");
    plotTemplate.SetDimensions(710, 710, true);
    plotTemplate.SetTransparent();
    plotTemplate[0].SetPosition(0., 0., 1., 1.);
    plotTemplate[0].SetFrameFill(10, 1001);
    plotTemplate[0].SetDefaultMarkerColors({kBlack, kBlue, kRed});
    plotTemplate[0].SetDefaultTextFont(43);
    plotTemplate[0].SetDefaultTextSize(24);
    plotTemplate[0].SetDefaultMarkerSize(1.);
    plotTemplate[0].SetDefaultLineWidth(3.);
    plotTemplate[0].SetDefaultDrawingOptionGraph(points);
    plotTemplate[0].SetTransparent();
    plotTemplate[0].SetMargins(0.05, 0.09, 0.09, 0.025);
    plotTemplate[0]['X'].SetTitleSize(28).SetTitleOffset(1.1).SetOppositeTicks().SetMaxDigits(3).SetNoExponent().SetTitleFont(63);
    plotTemplate[0]['Y'].SetTitleSize(28).SetTitleOffset(1.).SetOppositeTicks().SetMaxDigits(3).SetTitleFont(63);
    plotManager.AddPlotTemplate(plotTemplate);
  }  
  {
    Plot plotTemplate("lRaw", "PLOT_TEMPLATES");
    plotTemplate.SetDimensions(710, 710, true);
    plotTemplate.SetTransparent();
    plotTemplate[0].SetPosition(0., 0., 1., 1.);
    plotTemplate[0].SetFrameFill(10, 1001);
    plotTemplate[0].SetDefaultMarkerColors({kBlack, kBlue, kRed});
    plotTemplate[0].SetDefaultTextFont(43);
    plotTemplate[0].SetDefaultTextSize(24);
    plotTemplate[0].SetDefaultMarkerSize(1.);
    plotTemplate[0].SetDefaultLineWidth(3.);
    plotTemplate[0].SetDefaultDrawingOptionGraph(points);
    plotTemplate[0].SetTransparent();
    plotTemplate[0].SetMargins(0.05, 0.09, 0.16, 0.025);
    plotTemplate[0]['X'].SetTitleSize(28).SetTitleOffset(1.1).SetOppositeTicks().SetMaxDigits(3).SetNoExponent().SetTitleFont(63).SetLog(1);
    plotTemplate[0]['Y'].SetTitleSize(28).SetTitleOffset(1.5).SetOppositeTicks().SetMaxDigits(3).SetTitleFont(63).SetLog(1);
    plotManager.AddPlotTemplate(plotTemplate);
  }

  
}

//****************************************************************************************
/**
 * Defines the plots.
 */
//****************************************************************************************
void DefinePlots(PlotManager& plotManager)
{
//*****************************SINGLE HISTOS *******************************
  
  // -----------------------------------------------------------------------
  // -----------------******* INPUT SPECTRA  *********----------------------
  // -----------------------------------------------------------------------
  {
    Plot plot("genNewBin", "input", "myOwnLayout");      
    plot[1].AddData("hGenNB", "extrFileMC");        //[1] = pad!
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  } 


  // -----------------------------------------------------------------------
  // -----------------**********  SIGNAL  ************----------------------
  // -----------------------------------------------------------------------
  std::vector<string> ptLabels = {"0.1", "0.4", "0.5", "0.7", "0.8", "0.9", "1", "2", "3", "5", "18"};


  // vector<string> ptLabels = {"0.1", "0.3"};
  for(int iPt = 0; iPt < 10; iPt++){
    string ptLabel = ptLabels[iPt]+ " GeV/c < #it{p}_{T} < " + ptLabels[iPt+1] + " GeV/c";

    // ---Data----
    Plot plotData("signalData" + std::to_string(iPt), "signal", "lMyTemplate");   
    plotData[0]['X'].SetRange(0.35, 0.7);     //Eta-Rage
    // plotData[0]['X'].SetRange(0.05, 0.35);     //Pion-Rage
    plotData[0]['Y'].SetTitle("#it{N}");
    plotData[0]['X'].SetTitle("#it{m}_{ee #gamma}");
    plotData[1].AddData("hSignal" + std::to_string(iPt), "extrFileData", "Data");        //[1] = pad!
    plotData[1].AddData("fEtaPeak" + std::to_string(iPt), "extrFileData", "Gausian + pol2"); 
    plotData[1].AddData("hTrueProjection" + std::to_string(iPt), "extrFileMC", "MC true").SetColor(kOrange);
    plotData[1].AddText("ALICE work in progress // LHC24, LHC24f4d // #sqrt{s} = 13.6 TeV // #eta #rightarrow e^{+} e^{-} #gamma");
    plotData[1].AddLegend();
    plotData[1].AddText(ptLabel);
    plotManager.AddPlot(plotData);

    // ---reconstructed MC----
    Plot plotMC("signalMC" + std::to_string(iPt), "signal", "lMyTemplate");   
    plotMC[0]['X'].SetRange(0.35, 0.7);     //Eta-Rage
    // plotData[0]['X'].SetRange(0.05, 0.35);     //Pion-Rage 
    plotMC[0]['Y'].SetTitle("#it{N}");
    plotMC[0]['X'].SetTitle("#it{m}_{ee #gamma}");
    plotMC[1].AddData("hSignal"+ std::to_string(iPt), "extrFileMC", "rec. MC");        
    plotMC[1].AddData("fEtaPeak"+ std::to_string(iPt), "extrFileMC", "Gausian + pol2"); 
    plotMC[1].AddData("hTrueProjection"+ std::to_string(iPt), "extrFileMC", "MC true").SetColor(kOrange);
    plotMC[1].AddText("ALICE work in progress // LHC24f4d // #sqrt{s} = 13.6 TeV // #eta #rightarrow e^{+} e^{-} #gamma");
    plotMC[1].AddLegend();
    plotMC[1].AddText(ptLabel);
    plotManager.AddPlot(plotMC);
  }

 
  
    // string ptLabel = "
  // -----------------------------------------------------------------------
  // -----------------********  RAW SPECTRA **********----------------------
  // -----------------------------------------------------------------------
  {  
    Plot plot("rawSpecData", "rawSpec", "lMyTemplate");
    plot[1].AddData("hRawpT", "extrFileData", "R");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }  
  {  
    Plot plot("rawSpecMC", "rawSpec", "lMyTemplate");
    plot[1].AddData("hRawpT", "extrFileMC", "R");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }  
  {  
    Plot plot("trueRawSpec", "rawSpec", "lMyTemplate");
    plot[1].AddData("hTrueRawpT", "extrFileMC", "R");
    plot[1]['X'].SetLog();
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("genRawSpec", "rawSpec", "lMyTemplate");
    plot[1].AddData("hGenRawpT", "extrFileMC", "R");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("genDaliRawSpec", "rawSpec", "lMyTemplate");
    plot[1].AddData("hGenRawpTDalitz", "extrFileMC", "R");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }

  // -----------------------------------------------------------------------
  // ---------------------********** EFFI **********------------------------
  // -----------------------------------------------------------------------
  {  
    Plot plot("recEffi", "effi", "lMyTemplate");
    plot[1].AddData("hRecEffi", "effiFiles");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("trueEffi", "effi", "lMyTemplate");
    plot[1].AddData("hTrueEffi", "effiFiles");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("recDaliEffi", "effi", "lMyTemplate");
    plot[1].AddData("hRecEffiDalitz", "effiFiles");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("trueDaliEffi", "effi", "lMyTemplate");
    plot[1].AddData("hTrueEffiDalitz", "effiFiles");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }

  // -----------------------------------------------------------------------
  // -----------------*********  CORRECTED  **********----------------------
  // -----------------------------------------------------------------------
  {  
    Plot plot("recCorrData", "corrected", "lRaw");
    plot[1].AddData("hRecCorr", "corrFilesData", "hRecCorr");
    plot[1]['X'].SetLog();
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("recCorrMC", "corrected", "lRaw");
    plot[1].AddData("hRecCorr", "corrFilesMC", "hRecCorrMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("trueCorr", "corrected", "lRaw");
    plot[1].AddData("hTrueCorrNorm", "corrFilesMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("recDaliCorr", "corrected", "lRaw");
    plot[1].AddData("hRecCorrDalitz", "corrFilesMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("trueDaliCorr", "corrected", "lRaw");
    plot[1].AddData("hTrueCorrDalitz", "corrFilesMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }

  // -----------------------------------------------------------------------
  // -----------------*********  Checks **********----------------------
  // -----------------------------------------------------------------------
  {  
    Plot plot("significance", "checks", "lMyTemplate");
    plot[1].AddData("hSigni", "extrFileMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
  {  
    Plot plot("SoverB", "checks", "lMyTemplate");
    plot[1].AddData("hSB", "extrFileMC");
    // plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }

  // -----------------------------------------------------------------------
  //*****************************MULTI HISTOS *******************************
  // -----------------------------------------------------------------------

  {  
    Plot plot("trueVSgen", "closure", "lRaw");
    plot[1].AddData("hGenNB", "extrFileMC", " generated, rebinned").SetColor(kRed + 2).SetLegendID(2).SetMarkerSize(2).SetMarkerStyle(24);
    // plot[1].AddData("hTrueCorr", "corrFilesMC", "corrected true").SetColor(kBlue + 2);
    plot[1].AddData("hRecCorr", "corrFilesMC", "hRecCorrMC").SetColor(kBlue + 2);
    plot[1]['X'].SetLog();
    plot[1]['Y'].SetLog();
    plot[1].AddLegend();
    plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }

  {  
    Plot plot("corrSpectra", "closure", "lRaw");
    plot[1].AddData("hRecCorr", "corrFilesData", "hRecCorrData").SetColor(kBlue + 2);
    plot[1].AddData("hRecCorr", "corrFilesMC", "hRecCorrMC").SetColor(kRed + 2).SetLegendID(2);
    // plot[1]['X'].SetLog(0);
    plot[1]['Y'].SetLog();
    plot[1].AddLegend();
    plot[1].AddLegend();
    plotManager.AddPlot(plot);
  }
 
  // -----------------------------------------------------------------------
  // {  // -----------------------------------------------------------------------
  //   string plotName = "invMassDist";
  //   Plot plot(plotName, "examples", "lMyTemplate");
  //   plot[1].AddData("invMassDist", "pseudodata").SetScaleFactor(0.00001).SetColor(kBlue + 2);
  //   plotManager.AddPlot(plot);
  // }  // -----------------------------------------------------------------------
  // {  // -----------------------------------------------------------------------
  //   string plotName = "multDist";
  //   Plot plot(plotName, "examples", "lMyTemplate");
  //   plot[1].AddData("multDist", "pseudodata");
  //   plot[1].AddText("My Experiment // pp collisions, X TeV");
  //   plot[1]['Y'].SetLog();
  //   plotManager.AddPlot(plot);
  // }  // -----------------------------------------------------------------------
  // {  // -----------------------------------------------------------------------
  //   string plotName = "ptSpec";
  //   Plot plot(plotName, "examples", "lMyTemplate");
  //   plot[1].AddData("ptSpec", "pseudodata").SetColor(kRed + 2);
  //   plot[1]['X'].SetLog();
  //   plot[1]['Y'].SetLog();
  //   plotManager.AddPlot(plot);
  // }  // -----------------------------------------------------------------------
}
