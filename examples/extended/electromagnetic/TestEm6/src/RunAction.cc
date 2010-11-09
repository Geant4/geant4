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
// $Id: RunAction.cc,v 1.15 2010-11-09 21:30:47 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  for (G4int j=0; j<6; j++) histo[j] = 0;
  
#ifdef G4ANALYSIS_USE
 // Creating the analysis factory
 af = AIDA_createAnalysisFactory();
 
 if (af) { 
   // Creating the tree factory
   AIDA::ITreeFactory* tf = af->createTreeFactory();
 
   // Creating a tree mapped to an hbook file.
   G4bool readOnly  = false;
   G4bool createNew = true;
   G4String options = "";
   //tree = tf->create("testem6.hbook","hbook",readOnly,createNew,options);
   tree = tf->create("testem6.root", "root",readOnly,createNew,options);
   //tree = tf->create("testem6.XML" ,"XML"  ,readOnly,createNew,options);
   delete tf;
   
   if (tree) {   
     // Creating a histogram factory
     AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

     // Creating histograms
     histo[0] = hf->createHistogram1D("1","1/(1+(theta+[g]+)**2)",100, 0 ,1.);
     histo[1] = hf->createHistogram1D("2","log10(theta+ [g]+)",   100,-3.,1.);
     histo[2] = hf->createHistogram1D("3","log10(theta- [g]-)",   100,-3.,1.);
     histo[3] = hf->createHistogram1D("4","log10(theta+ [g]+ -theta- [g]-)",
                                                                  100,-3.,1.);
     histo[4] = hf->createHistogram1D("5","xPlus" ,100,0.,1.);
     histo[5] = hf->createHistogram1D("6","xMinus",100,0.,1.);
  
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
  ProcCounter = new ProcessesCount;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
  //does the process  already encounted ?
  size_t nbProc = ProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));
  
  (*ProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
  //total number of process calls
  G4double countTot = 0.;
  G4cout << "\n Number of process calls --->";
  for (size_t i=0; i< ProcCounter->size();i++) {
	G4String procName = (*ProcCounter)[i]->GetName();
	if (procName != "Transportation") {
	  G4int count    = (*ProcCounter)[i]->GetCounter(); 
	  G4cout << "\t" << procName << " : " << count;
	  countTot += count;
	}
  }
  G4cout << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
