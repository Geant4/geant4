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
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALRunAction.cc,v 1.3 2002-12-17 15:45:42 pmendez Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE

#include "FCALRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "Randomize.hh"

#include <AIDA/AIDA.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALRunAction::FCALRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALRunAction::~FCALRunAction()
{

  cleanHisto();

}


void FCALRunAction::bookHisto(){

  AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
  AIDA::ITreeFactory* tf = af->createTreeFactory();

 // Creating a tree mapped to an hbook file.
  tree = tf->create("fcal.paw", "hbook", false, true);

 // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

 // Creating a ntuple factory.

  AIDA::ITupleFactory* nf = af->createTupleFactory(*tree);

 // Creating ntuples

  ntuple[1] = nf->create("100","Number Out of World","float OutOfWorld, i, j");
   
  ntuple[2] = nf->create("200","Secondary Info","float Secondary, i, j");
  
  ntuple[3] = nf->create("300","Energy Deposits","float EmEdep, HadEdep");

  // Creating Histograms
  
  histo[1] = hf->createHistogram1D("1","Number of OutOfWorld",100, 0., 100.);
  histo[2] = hf->createHistogram1D("2","Number of Secondaries",100, 0., 100.);
  histo[3] = hf->createHistogram1D("3","Electromagnetic Energy / MeV",100, 0., 100.);
  histo[4] = hf->createHistogram1D("4","Hadronic Energy / MeV",100, 0., 100.);

  delete hf;
  delete tf;
  delete af;     
  delete nf;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCALRunAction::cleanHisto()
{

  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
  
  delete tree;
 

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if (aRun->GetRunID() == 0) bookHisto();

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALRunAction::EndOfRunAction(const G4Run* )
{
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
}

#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






