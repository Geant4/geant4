//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
fileName("exgps.aida"),fileType("xml"),analysisFactory(0), hFactory(0), tFactory(0),
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
  if (analysisFactory)
  {
    IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
    if (pf) return pf->create("Plotter");
  }
  return 0;
}


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
  enerHisto->fill(e/MeV,w);
  posiXY->fill(x/cm,y/cm,w);
  posiXZ->fill(x/cm,z/cm,w);
  posiYZ->fill(y/cm,z/cm,w);
  anglCTP->fill(p/deg,std::cos(t),w);
  anglTP->fill(p/deg,t/deg,w);

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

  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory){
    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create(fileName,fileType,false,true,"compress=yes");
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
  }

  // Create a Tuple

  if (tFactory)
  {
     tuple = tFactory->create("MyTuple","MyTuple","string Pname, double Energy, X, Y, Z, Theta, Phi, Weight","");
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
    delete tFactory;
    delete hFactory;
    //    G4cout << "Warning: Geant4 will NOT continue unless you close the JAS-AIDA window." << G4endl;
    delete analysisFactory;
  }
  //  dispose();
}


#endif // G4ANALYSIS_USE










