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
/// \file eventgenerator/exgps/src/exGPSAnalysisManager.cc
/// \brief Implementation of the exGPSAnalysisManager class
//
#ifdef G4ANALYSIS_USE
#include <AIDA/AIDA.h>
#include "exGPSAnalysisManager.hh"
#include "exGPSAnalysisMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

exGPSAnalysisManager* exGPSAnalysisManager::fInstance = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::exGPSAnalysisManager(): 
fIleName("exgps.aida"),fIleType("xml"),fAnalysisFactory(0),
fTree(0),fPlotter(0),fMinpos(-10.),fMaxpos(10),fMineng(0.),fMaxeng(1000.),
fEnerHisto(0),fPosiXY(0),fPosiXZ(0),fPosiYZ(0),fAnglCTP(0),fAnglTP(0),fTuple(0)
{
  // Define the messenger and the analysis system
  fAnalysisMessenger = new exGPSAnalysisMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager::~exGPSAnalysisManager() {
  delete fAnalysisMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exGPSAnalysisManager* exGPSAnalysisManager::GetInstance ()
{
  if (fInstance == 0) fInstance = new exGPSAnalysisManager();
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisManager::Dispose()
{
  if (fInstance != 0)
  {
    delete fInstance;
    fInstance = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exGPSAnalysisManager::Fill(G4String pname, G4double e, 
                                G4double x, G4double y, G4double z,
                                G4double t, G4double p, G4double w)
{
  if(fEnerHisto) {
    fEnerHisto->fill(e/MeV,w);
    fPosiXY->fill(x/cm,y/cm,w);
    fPosiXZ->fill(x/cm,z/cm,w);
    fPosiYZ->fill(y/cm,z/cm,w);
    fAnglCTP->fill(p/deg,std::cos(t),w);
    fAnglTP->fill(p/deg,t/deg,w);
  }

  if (fPlotter) fPlotter->refresh();

  // Fill the fTuple

  if (fTuple) {
    fTuple->fill(0,pname);
    fTuple->fill(1,e/MeV);
    fTuple->fill(2,x/cm);
    fTuple->fill(3,y/cm);
    fTuple->fill(4,z/cm);
    fTuple->fill(5,t/deg);
    fTuple->fill(6,p/deg);
    fTuple->fill(7,w);
    fTuple->addRow();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the initialization so that the histograms have
   always the right dimensions depending from the detector geometry
*/
void exGPSAnalysisManager::BeginOfRun() 
{ 
  fTree = 0;
  fPlotter = 0;
  fEnerHisto =0;
  fPosiXY = fPosiXZ = fPosiYZ = fAnglCTP =fAnglTP = 0;
  fTuple = 0;

  //Hooking an AIDA compliant analysis system.
  fAnalysisFactory = AIDA_createAnalysisFactory();
  if(!fAnalysisFactory) //Have to check that, it can fail.
  {
    G4cout << "exGPSAnalysisManager::BeginOfRun: can't get AIDA." << G4endl; 
    return;
  }
  AIDA::ITreeFactory* treeFactory = fAnalysisFactory->createTreeFactory();
  if(treeFactory) 
  {
    fTree = treeFactory->create(fIleName,fIleType,false,true,"compress=yes");
    if(!fTree) //Have to check that, it can fail.
    {
      G4cout << "exGPSAnalysisManager::BeginOfRun:"
             << " can't create the AIDA::ITree : " << fIleName << G4endl; 
      return;
    }

    delete treeFactory; // Will not delete the ITree.
  }

  AIDA::IHistogramFactory* hFactory =
    fAnalysisFactory->createHistogramFactory(*fTree);
  if (hFactory)
  {
    // Create the energy histogram
    fEnerHisto = hFactory->createHistogram1D("Source Energy Spectrum",100,
                                                                                                                      fMineng,fMaxeng);

    // Create some 2d histos 
    fPosiXY = hFactory->createHistogram2D("Source X-Y distribution",100,
                              fMinpos/cm,fMaxpos/cm,100,fMinpos/cm,fMaxpos/cm);
    fPosiXZ = hFactory->createHistogram2D("Source X-Z distribution",100,
                              fMinpos/cm,fMaxpos/cm,100,fMinpos/cm,fMaxpos/cm);
    fPosiYZ = hFactory->createHistogram2D("Source Y-Z distribution",100,
                              fMinpos/cm,fMaxpos/cm,100,fMinpos/cm,fMaxpos/cm);
    fAnglCTP = hFactory->createHistogram2D("Source phi-std::cos(theta) distribution",
                              360,0,360,100, -1, 1);
    fAnglTP = hFactory->createHistogram2D("Source phi-theta distribution",
                              360,0,360,180,0,180);
    delete hFactory;
  }

  // Create a Tuple

  AIDA::ITupleFactory* tFactory = fAnalysisFactory->createTupleFactory(*fTree);
  if (tFactory)
  {
     fTuple = tFactory->create("MyTuple","MyTuple",
                     "string Pname, double Energy, X, Y, Z, Theta, Phi, Weight","");
     delete tFactory;
  }
  
  AIDA::IPlotterFactory* pf = fAnalysisFactory->createPlotterFactory(0,0);
  if(pf) 
  {
    fPlotter = pf->create();
    if (fPlotter)
    {
       fPlotter->createRegions(2,3);
       if(fEnerHisto) fPlotter->region(0)->plot(*fEnerHisto);
       if(fPosiXY) fPlotter->region(1)->plot(*fPosiXY);
       if(fPosiXZ) fPlotter->region(2)->plot(*fPosiXZ);
       if(fPosiYZ) fPlotter->region(3)->plot(*fPosiYZ);
       if(fAnglCTP) fPlotter->region(4)->plot(*fAnglCTP);
       if(fAnglTP) fPlotter->region(5)->plot(*fAnglTP);
       fPlotter->show();
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
  if (fAnalysisFactory)
  {
    if (!fTree->commit()) G4cout << "Commit failed: no AIDA file produced!"
                                                                                                                                    << G4endl;
    delete fTree;
    fTree = 0;

    delete fPlotter;
    fPlotter = 0;

    delete fAnalysisFactory; //cleanup all AIDA related things.
    fAnalysisFactory = 0;

    fEnerHisto = 0;
    fPosiXY = fPosiXZ = fPosiYZ = fAnglCTP =fAnglTP = 0;
    fTuple = 0;

  }
}


#endif // G4ANALYSIS_USE
