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
// $Id: RunAction.cc,v 1.9 2005/12/06 11:44:55 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
:af(0), tree(0)
{
  histo[0] = 0;
  
#ifdef G4ANALYSIS_USE 
 // Creating the analysis factory
 af = AIDA_createAnalysisFactory();
 
 if (af) {
   // Creating the tree factory
   AIDA::ITreeFactory* tf = af->createTreeFactory();
 
   // Creating a tree mapped to an hbook file.
   G4bool readOnly  = false;
   G4bool createNew = true;
   G4String options =  "--noErrors uncompress";
   tree = tf->create("testem4.hbook","hbook",readOnly,createNew, options);
   //tree = tf->create("testem4.root","root",readOnly,createNew, options);
   //tree = tf->create("testem4.XML" ,"XML" ,readOnly,createNew, options);
   delete tf;
   
   if (tree) {
     // Creating a histogram factory
     AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

     // Creating the histogram
     histo[0]=hf->createHistogram1D
                         ("1","total energy deposit in C6F6(MeV)",100,0.,10.);

     delete hf;
     G4cout << "\n----> Histogram tree is opened" << G4endl;
   }
 }
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
#ifdef G4ANALYSIS_USE
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
    
  delete tree;
  delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
