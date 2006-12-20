//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifdef G4ANALYSIS_USE

#ifdef G4ANALYSIS_USE_AIDA
#include <AIDA/AIDA.h>
#endif

#ifdef G4ANALYSIS_USE_ROOT
#include "TROOT.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#endif


#include "exGPSAnalysisManager.hh"
#include "exGPSAnalysisMessenger.hh"

#include "G4UnitsTable.hh"

exGPSAnalysisManager* exGPSAnalysisManager::instance = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::exGPSAnalysisManager(): fileName("exgps"),
#ifdef G4ANALYSIS_USE_AIDA
fileType("xml"),analysisFactory(0), hFactory(0), tFactory(0),
#endif
minpos(-10.),maxpos(10),mineng(0.),maxeng(1000.)
{
  // Define the messenger and the analysis system
  analysisMessenger = new exGPSAnalysisMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::~exGPSAnalysisManager() {

  delete analysisMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE_AIDA

IHistogramFactory* exGPSAnalysisManager::getHistogramFactory()
{
  return hFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


ITupleFactory* exGPSAnalysisManager::getTupleFactory()
{
  return tFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

IPlotter* exGPSAnalysisManager::createPlotter()
{
#ifdef JAIDA_HOME 
  if (analysisFactory)
  {
    IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
    if (pf) return pf->create("Plotter");
  }
#endif
  return 0;
}

#endif
////////////////////////////////////////////////////////////////////////////////
//
exGPSAnalysisManager* exGPSAnalysisManager::getInstance ()
{
  if (instance == 0) instance = new exGPSAnalysisManager();
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisManager::dispose()
{
  if (instance != 0)
  {
    delete instance;
    instance = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisManager::Fill(G4double pid, G4double e, 
				G4double x, G4double y, G4double z,
				G4double t, G4double p, G4double w)
{
#ifdef G4ANALYSIS_USE_AIDA
  enerHisto->fill(e/MeV,w);
  posiXY->fill(x/cm,y/cm,w);
  posiXZ->fill(x/cm,z/cm,w);
  posiYZ->fill(y/cm,z/cm,w);
  anglCTP->fill(p/deg,std::cos(t),w);
  anglTP->fill(p/deg,t/deg,w);

  if (plotter) plotter->refresh();

  // Fill the tuple

  if (tuple) {
    tuple->fill(0,pid);
    tuple->fill(1,e/MeV);
    tuple->fill(2,x/cm);
    tuple->fill(3,y/cm);
    tuple->fill(4,z/cm);
    tuple->fill(5,t/deg);
    tuple->fill(6,p/deg);
    tuple->fill(7,w);
    
    tuple->addRow();    
  }
#endif

#ifdef G4ANALYSIS_USE_ROOT
  enerHistoroot->Fill(e/MeV,w);
  posiXYroot->Fill(x/cm,y/cm,w);
  posiXZroot->Fill(x/cm,z/cm,w);
  posiYZroot->Fill(y/cm,z/cm,w);
  anglCTProot->Fill(p/deg,std::cos(t),w);
  anglTProot->Fill(p/deg,t/deg,w);

  // Fill the tuple

  if (tupleroot) {
    tupleroot->Fill(pid,e/MeV,x/cm,y/cm,z/cm,t/deg,p/deg,w);
  }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void exGPSAnalysisManager::BeginOfRun() 
{ 
#ifdef G4ANALYSIS_USE_AIDA
  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory){
    G4String AIDAfileName         = fileName + G4String(".aida");
    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create(AIDAfileName,fileType,false,true,"compress=yes");
    hFactory = analysisFactory->createHistogramFactory(*tree);
    tFactory = analysisFactory->createTupleFactory(*tree);
    delete treeFactory; // Will not delete the ITree.
  }
  //
  enerHisto =0;
  posiXY = posiXZ = posiYZ = anglCTP =anglTP = 0;
  plotter = 0;
  tuple = 0;
  //
  if (hFactory)
  {
    // Create the energy histogram
    enerHisto = hFactory->createHistogram1D("Source Energy Spectrum",100,mineng,maxeng);

    // Create some 2d histos 
    posiXY = hFactory->createHistogram2D("Source X-Y distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
    posiXZ = hFactory->createHistogram2D("Source X-Z distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
    posiYZ = hFactory->createHistogram2D("Source Y-Z distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
    anglCTP = hFactory->createHistogram2D("Source phi-std::cos(theta) distribution",360,0,360
						   , 100, -1, 1);
    anglTP = hFactory->createHistogram2D("Source phi-theta distribution",360,0,360
						  ,180,0,180);
#ifdef JAIDA_HOME
    plotter = createPlotter();

    if (plotter)
    {
       plotter->createRegions(2,3);
       plotter->region(0)->plot(*enerHisto);
       plotter->region(1)->plot(*posiXY);
       plotter->region(2)->plot(*posiXZ);
       plotter->region(3)->plot(*posiYZ);
       plotter->region(4)->plot(*anglCTP);
       plotter->region(5)->plot(*anglTP);
       plotter->show();
     }
#endif

  }

  // Create a Tuple

  if (tFactory)
  {
     tuple = tFactory->create("10","10","std::double Pid, Energy, X, Y, Z, Theta, Phi, Weight","");
  }
#endif

#ifdef G4ANALYSIS_USE_ROOT
  new TApplication("App", ((int *)0), ((char **)0));
  G4String ROOTfileName = fileName + G4String(".root");
  hfileroot = new TFile(ROOTfileName.c_str() ,"RECREATE","ROOT file for exGPS");
  
  // Create the energy histogram
  enerHistoroot =new TH1D("h1000","Source Energy Spectrum",100,mineng,maxeng);

  // Create some 2d histos 
  posiXYroot = new TH2D("h2001","Source X-Y distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
  posiXZroot = new TH2D("h2002","Source X-Z distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
  posiYZroot = new TH2D("h2003","Source Y-Z distribution",100,minpos/cm,maxpos/cm
					    ,100,minpos/cm,maxpos/cm);
  anglCTProot =new TH2D("h2004","Source phi-std::cos(theta) distribution",360,0,360
						   , 100, -1, 1);
  anglTProot = new TH2D("h2005","Source phi-theta distribution",360,0,360
						  ,180,0,180);
  
  tupleroot =  new TNtuple("t10","Event by event data","Pname:X:Y:Z:Theta:Phi:Weight"); 
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void exGPSAnalysisManager::EndOfRun() 
{
#ifdef G4ANALYSIS_USE_AIDA
  if (analysisFactory)
  {
    if (!tree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
    delete tree;
    delete tFactory;
    delete hFactory;
    //    G4cout << "Warning: Geant4 will NOT continue unless you close the JAS-AIDA window." << G4endl;
    delete analysisFactory;
  }
  //  dispose();
#endif

#ifdef G4ANALYSIS_USE_ROOT
  G4cout << "ROOT: files writing..." << G4endl;
  hfileroot->Write();
  G4cout << "ROOT: files closing..." << G4endl;
  hfileroot->Close();
#endif
}

#endif // G4ANALYSIS_USE










