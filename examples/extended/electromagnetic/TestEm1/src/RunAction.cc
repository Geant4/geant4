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
// $Id: RunAction.cc,v 1.3 2003/11/07 15:38:28 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 08.03.01 Hisaya: adapted for STL   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#ifdef G4ANALYSIS_USE
 #include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
  : ProcCounter(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
 cleanHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
#ifdef G4ANALYSIS_USE
 // Creating the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 AIDA::ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create("testem1.paw", "hbook", readOnly, createNew);

 // Creating a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);


 // booking histograms
 histo[0] = hf->createHistogram1D("1","track length (mm) of a charged particle",
                         100,0.,50*cm);
 histo[1] = hf->createHistogram1D("2","Nb of steps per track (charged particle)",
                         100,0.,100.);
 histo[2] = hf->createHistogram1D("3","step length (mm) charged particle",
                         100,0.,10*mm);
		       
 delete hf;
 delete tf;
 delete af;		       
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::cleanHisto()
{
#ifdef G4ANALYSIS_USE
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
 
  delete tree;
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  NbOfTraks0 = 0; NbOfTraks1 = 0; NbOfSteps0 = 0; NbOfSteps1 = 0;
  edep = 0.0;
  ProcCounter = new ProcessesCount;
     
  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
    
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
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

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents)
    { //nb of tracks and steps per event
      G4double dNbOfEvents = double(NbOfEvents);
    
      std::ios::fmtflags mode = G4cout.flags();
      G4cout.setf(std::ios::fixed,std::ios::floatfield);

      G4int  prec = G4cout.precision(4);
      
      G4cout << "\n nb tracks/event"
             << "   neutral: " << std::setw(10) << NbOfTraks0/dNbOfEvents
             << "   charged: " << std::setw(10) << NbOfTraks1/dNbOfEvents
             << "\n nb  steps/event"
             << "   neutral: " << std::setw(10) << NbOfSteps0/dNbOfEvents
             << "   charged: " << std::setw(10) << NbOfSteps1/dNbOfEvents
             << G4endl;
	     
      G4cout << "\n total energy deposit: " 
             << G4BestUnit(edep/dNbOfEvents, "Energy") << G4endl;
      
      //frequency of processes call       
      G4cout << "\n nb of process calls per event: \n   ";       
      for (size_t i=0; i< ProcCounter->size();i++)
           G4cout << std::setw(12) << (*ProcCounter)[i]->GetName();
           
      G4cout << "\n   ";       
      for (size_t j=0; j< ProcCounter->size();j++)
      G4cout << std::setw(12) << ((*ProcCounter)[j]->GetCounter())
                                                               /dNbOfEvents;
      G4cout << G4endl;    
                         
      G4cout.setf(mode,std::ios::floatfield);
      G4cout.precision(prec);       
    }         

  // delete and remove all contents in ProcCounter 
  while (ProcCounter->size()>0){
    OneProcessCount* aProcCount=ProcCounter->back();
    ProcCounter->pop_back();
    delete aProcCount;
  }
  delete ProcCounter;
                             
  //draw the events
  if (G4VVisManager::GetConcreteInstance()) 
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
