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
// $Id: Em1RunAction.cc,v 1.18 2002-12-11 14:19:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 08.03.01 Hisaya: adapted for STL   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em1RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "Randomize.hh"

#ifndef G4NOHIST
 #include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1RunAction::Em1RunAction()
  : ProcCounter(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1RunAction::~Em1RunAction()
{
 cleanHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1RunAction::bookHisto()
{
#ifndef G4NOHIST
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

void Em1RunAction::cleanHisto()
{
#ifndef G4NOHIST
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
 
  delete tree;
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  NbOfTraks0 = 0; NbOfTraks1 = 0; NbOfSteps0 = 0; NbOfSteps1 = 0;
  ProcCounter = new ProcessesCount;
     
  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
    
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1RunAction::CountProcesses(G4String procName)
{
   //does the process  already encounted ?
   size_t nbProc = ProcCounter->size();
   size_t i = 0;
   while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
   if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));

   (*ProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents)
    { //nb of tracks and steps per event
      G4double dNbOfEvents = double(NbOfEvents);
    
#ifdef G4USE_STD_NAMESPACE
      G4std::ios::fmtflags mode = G4cout.flags();
      G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#else 
      G4long mode = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#endif

      G4int  prec = G4cout.precision(4);
      
      G4cout << "\n nb tracks/event"
             << "   neutral: " << G4std::setw(10) << NbOfTraks0/dNbOfEvents
             << "   charged: " << G4std::setw(10) << NbOfTraks1/dNbOfEvents
             << "\n nb  steps/event"
             << "   neutral: " << G4std::setw(10) << NbOfSteps0/dNbOfEvents
             << "   charged: " << G4std::setw(10) << NbOfSteps1/dNbOfEvents
             << G4endl;
      
      //frequency of processes call       
      G4cout << "\n nb of process calls per event: \n   ";       
      for (size_t i=0; i< ProcCounter->size();i++)
           G4cout << G4std::setw(12) << (*ProcCounter)[i]->GetName();
           
      G4cout << "\n   ";       
      for (size_t j=0; j< ProcCounter->size();j++)
      G4cout << G4std::setw(12) << ((*ProcCounter)[j]->GetCounter())
                                                               /dNbOfEvents;
      G4cout << G4endl;    
                         
      G4cout.setf(mode,G4std::ios::floatfield);
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
