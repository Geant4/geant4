// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1RunAction.cc,v 1.9 2001-03-08 14:57:43 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// 08.03.01 Hisaya: adapted for STL   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1RunAction.hh"
#include "Em1RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "Randomize.hh"

#ifndef G4NOHIST
 #include "CLHEP/Hist/HBookFile.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1RunAction::Em1RunAction()
  : ProcCounter(0), saveRndm (1),
    runMessenger(new Em1RunActionMessenger(this))
{
#ifndef G4NOHIST
  hbookManager = NULL;
#endif 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1RunAction::~Em1RunAction()
{
 cleanHisto();
 delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunAction::bookHisto()
{
#ifndef G4NOHIST
  hbookManager = new HBookFile("testem1.histo", 68);

  // booking histograms
  histo[0] = hbookManager->histogram("track length (mm) of a charged particle",100,0.,50*cm);
  histo[1] = hbookManager->histogram("Nb of steps per track (charged particle)",100,0.,100.);
  histo[2] = hbookManager->histogram("step length (mm) charged particle",100,0.,10*mm);
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunAction::cleanHisto()
{
#ifndef G4NOHIST
  delete [] histo;
  delete hbookManager;
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
    }  
  
  NbOfTraks0 = 0; NbOfTraks1 = 0; NbOfSteps0 = 0; NbOfSteps1 = 0;
  ProcCounter = new ProcessesCount;
     
  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
    
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunAction::CountProcesses(G4String procName)
{
   //does the process  already encounted ?
   size_t nbProc = ProcCounter->size();
   size_t i = 0;
   while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
   if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));

   (*ProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents)
    { //nb of tracks and steps per event
      G4double dNbOfEvents = double(NbOfEvents);
    
      G4long oldform = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
      G4int  oldprec = G4cout.precision(2);
      
      G4cout << "\n nb tracks/event"
                      << "   neutral: " << G4std::setw(7) << NbOfTraks0/dNbOfEvents
                      << "   charged: " << G4std::setw(7) << NbOfTraks1/dNbOfEvents
             << "\n nb  steps/event"
                      << "   neutral: " << G4std::setw(7) << NbOfSteps0/dNbOfEvents
                      << "   charged: " << G4std::setw(7) << NbOfSteps1/dNbOfEvents
             << G4endl;
      
      //frequency of processes call       
      G4cout << "\n nb of process calls per event: \n   ";       
      for (G4int i=0; i< ProcCounter->size();i++)
           G4cout << G4std::setw(9) << (*ProcCounter)[i]->GetName();
           
      G4cout << "\n   ";       
      for (G4int j=0; j< ProcCounter->size();j++)
           G4cout << G4std::setw(9) << ((*ProcCounter)[j]->GetCounter())/dNbOfEvents;
      G4cout << G4endl;    
                         
      G4cout.setf(oldform,G4std::ios::floatfield);
      G4cout.precision(oldprec);       
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

  // save Rndm status
  if (saveRndm > 0)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("endOfRun.rndm");
    }
    
#ifndef G4NOHIST     
  // writing histogram file
  hbookManager->write();
#endif               
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
