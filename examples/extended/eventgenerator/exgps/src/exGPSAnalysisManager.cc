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

#include <AIDA/AIDA.h>

#include "exGPSAnalysisManager.hh"
#include "exGPSAnalysisMessenger.hh"

#include "G4UnitsTable.hh"

exGPSAnalysisManager* exGPSAnalysisManager::instance = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::exGPSAnalysisManager(): 
fileName("exgps.aida"),fileType("xml"),analysisFactory(0), tree(0),plotter(0),
minpos(-10.),maxpos(10),mineng(0.),maxeng(1000.),
enerHisto(0),posiXY(0),posiXZ(0),posiYZ(0),anglCTP(0),anglTP(0),tuple(0)
{
  // Define the messenger and the analysis system
  analysisMessenger = new exGPSAnalysisMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::~exGPSAnalysisManager() {

  delete analysisMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

void exGPSAnalysisManager::Fill(G4String pname, G4double e, 
				G4double x, G4double y, G4double z,
				G4double t, G4double p, G4double w)
{
  if(enerHisto) {
    enerHisto->fill(e/MeV,w);
    posiXY->fill(x/cm,y/cm,w);
    posiXZ->fill(x/cm,z/cm,w);
    posiYZ->fill(y/cm,z/cm,w);
    anglCTP->fill(p/deg,std::cos(t),w);
    anglTP->fill(p/deg,t/deg,w);
  }

  if (plotter) plotter->refresh();

  // Fill the tuple

  if (tuple) {
    tuple->fill(0,pname);
    tuple->fill(1,e/MeV);
    tuple->fill(2,x/cm);
    tuple->fill(3,y/cm);
    tuple->fill(4,z/cm);
    tuple->fill(5,t/deg);
    tuple->fill(6,p/deg);
    tuple->fill(7,w);
    
    tuple->addRow();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void exGPSAnalysisManager::BeginOfRun() 
{ 
  tree = 0;
  plotter = 0;

  enerHisto =0;
  posiXY = posiXZ = posiYZ = anglCTP =anglTP = 0;
  tuple = 0;

  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(!analysisFactory) //Have to check that, it can fail.
  {
    G4cout << "exGPSAnalysisManager::BeginOfRun: can't get AIDA." << G4endl; 
    return;
  }
  AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
  if(treeFactory) 
  {
    tree = treeFactory->create(fileName,fileType,false,true,"compress=yes");
    if(!tree) //Have to check that, it can fail.
    {
      G4cout << "exGPSAnalysisManager::BeginOfRun:"
             << " can't create the AIDA::ITree : " << fileName << G4endl; 
      return;
    }

    delete treeFactory; // Will not delete the ITree.
  }

  AIDA::IHistogramFactory* hFactory = analysisFactory->createHistogramFactory(*tree);
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
    delete hFactory;
  }

  // Create a Tuple

  AIDA::ITupleFactory* tFactory = analysisFactory->createTupleFactory(*tree);
  if (tFactory)
  {
     tuple = tFactory->create("MyTuple","MyTuple","string Pname, double Energy, X, Y, Z, Theta, Phi, Weight","");
     delete tFactory;
  }
  
  AIDA::IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
  if(pf) 
  {
    plotter = pf->create();
    if (plotter)
    {
       plotter->createRegions(2,3);
       if(enerHisto) plotter->region(0)->plot(*enerHisto);
       if(posiXY) plotter->region(1)->plot(*posiXY);
       if(posiXZ) plotter->region(2)->plot(*posiXZ);
       if(posiYZ) plotter->region(3)->plot(*posiYZ);
       if(anglCTP) plotter->region(4)->plot(*anglCTP);
       if(anglTP) plotter->region(5)->plot(*anglTP);
       plotter->show();
    }
    delete pf;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void exGPSAnalysisManager::EndOfRun() 
{
  if (analysisFactory)
  {
    if (!tree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
    delete tree;
    tree = 0;
    //    G4cout << "Warning: Geant4 will NOT continue unless you close the JAS-AIDA window." << G4endl;

    delete plotter;
    plotter = 0;

    delete analysisFactory; //cleanup all AIDA related things.
    analysisFactory = 0;

    enerHisto = 0;
    posiXY = posiXZ = posiYZ = anglCTP =anglTP = 0;
    tuple = 0;

  }
}


#endif // G4ANALYSIS_USE










