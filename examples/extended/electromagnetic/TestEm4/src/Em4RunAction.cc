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
//
// $Id: Em4RunAction.cc,v 1.12 2002-05-29 15:32:28 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em4RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "Randomize.hh"

#ifndef G4NOHIST
 #include "AIDA/IAnalysisFactory.h"
 #include "AIDA/ITreeFactory.h"
 #include "AIDA/ITree.h"
 #include "AIDA/IHistogramFactory.h"
 #include "AIDA/IHistogram1D.h"
 #include "AIDA/IAxis.h"
 #include "AIDA/IAnnotation.h"
 #include "AIDA/ITupleFactory.h"
 #include "AIDA/ITuple.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em4RunAction::Em4RunAction()
{
  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em4RunAction::~Em4RunAction()
{
#ifndef G4NOHIST
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)

  delete tree;
  delete [] histo;
#endif  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em4RunAction::bookHisto()
{
#ifndef G4NOHIST 
 // Creating the analysis factory
 IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 tree = tf->create("testem4.paw", false, false, "hbook");

 // Creating a histogram factory, whose histograms will be handled by the tree
 IHistogramFactory* hf = af->createHistogramFactory(*tree);
 
 // Creating the histogram
 histo[0]=hf->create1D("eDep","total energy deposit in C6F6 (MeV)",100,0.,10.);
 
 delete hf;
 delete tf;
 delete af;
#endif      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em4RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em4RunAction::EndOfRunAction(const G4Run* aRun)
{  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  // show Rndm status
  HepRandom::showEngineStatus();         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
