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
// $Id: FCALAnalysisManager.cc
// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
// 12 Feb 2003 Patricia Mendez      Created
// -------------------------------------------------------------------
#include <stdlib.h>
#include "G4VProcess.hh"
#include <fstream>
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#include "FCALAnalysisManager.hh"
#include "FCALAnalysisMessenger.hh"


#include "G4Step.hh"

FCALAnalysisManager* FCALAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager::FCALAnalysisManager()
  :outputFileName("fcal.his"),analysisFactory(0), tree(0),histogramFactory(0) 
   //  OutOfWorld(0), Secondary(0), EmEdep(0), HadEdep(0)  
{

  analisysMessenger = new FCALAnalysisMessenger(this);

  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      // tree = treeFactory->create(); // Tree in memory.
      //porebbe essere qua il memory leak//

      tree = treeFactory->create(outputFileName,"hbook",false,true);

      delete treeFactory; // Will not delete the ITree.
      histogramFactory = analysisFactory->createHistogramFactory(*tree);  
      tupleFactory = analysisFactory->createTupleFactory(*tree);
    }
  }

  G4cout << "FCALAnalysisManager created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager::~FCALAnalysisManager() 
{
  delete histogramFactory;
  histogramFactory=0;

  delete analysisFactory;
  analysisFactory = 0;

  delete tupleFactory;
  tupleFactory=0;

  delete instance;

  G4cout << "FCALAnalysisManager delete" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager* FCALAnalysisManager::getInstance()

{
  if (instance == 0) {instance = new FCALAnalysisManager;}
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::book()
{
    // Book histograms

  histo_1 = histogramFactory->createHistogram1D("1","Number of Out Of World", 100,0.,10.); 
  histo_2 = histogramFactory->createHistogram1D("2","Number of Secondaries", 100,0.,100.);
  histo_3 = histogramFactory->createHistogram1D("3","Electromagnetic Energy/MeV", 100,0.,100.);
  histo_4 = histogramFactory->createHistogram1D("4","hadronic Energy/MeV", 100,10.,60.);
 
  // Create a tuple 

  //  tuple = tupleFactory->create("FCAL","FCAL","energy counts");

  // ntuple_1 = tupleFactory->create("100","Number Out of World","float OutOfWorld, i, j");
   
  // ntuple_2 = tupleFactory->create("200","Secondary Info","float Secondary, i, j");
  
  // ntuple_3 = tupleFactory->create("300","Energy Deposits","float EmEdep, HadEdep");
} 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::finish()
{

  if(tree) {
    tree->commit(); // Write histos and tuple in file. 
    tree->close();
  }

}

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
//void FCALAnalysisManager::analyseEnergyDep(G4double Edep)
//void FCALAnalysisManager::analyseEnergyDep(G4double HadEdep)
//{
//
//
//  histo_1->fill(OutOfWorld);
//  histo_2->fill(Secondary);
//  histo_3->fill(EmEdep);
//  histo_4->fill(HadEdep);
//
// 
//
//  
//  if(ntuple_1) {
//    ntuple_1->fill(0,OutOfWorld);
//    ntuple_1->fill(1,1);
//    ntuple_1->addRow();
//  }
//
//  if(ntuple_2) {
//    ntuple_2->fill(0,Secondary);
//    ntuple_2->fill(1,1);
//    ntuple_2->addRow();
//  }
//
//  if(ntuple_3) {
//    ntuple_3->fill(0,EmEdep);
//    ntuple_3->fill(1,1);
//    ntuple_3->addRow();
//  }
//
//  
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::SetOutputFileName(G4String newName)
{

  outputFileName = newName;

}

#endif










