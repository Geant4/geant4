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
//
//
#include "G4tbbWorkerRunManager.hh"

G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType 
 G4tbbWorkerRunManager::instancesList=
   G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType();

G4tbbWorkerRunManager::G4tbbWorkerRunManager() 
   : fWorkerContext(0)
{
  // Register this RM - needed to eventually free/destroy resources
  instancesList.push(this);

  //  InitSlave();
  fWorkerContext= new G4WorkerThread(); 

  fWorkerContext->BuildGeometryAndPhysicsVector(); 
}

G4tbbWorkerRunManager::~G4tbbWorkerRunManager() 
{
  fWorkerContext->DestroyGeometryAndPhysicsVector(); 

  //  These would do something like:  
  //    tbbSlaveDestroyGeometryAndPhysicsVector(); 

  delete fWorkerContext; 
  fWorkerContext= 0;
}

void G4tbbWorkerRunManager::ProcessOneEvent( G4int i_event)
{
   // In MT the seeds are set first

   // Do the rest
   currentEvent = GenerateEvent(i_event);
   eventManager->ProcessOneEvent(currentEvent);
   AnalyzeEvent(currentEvent);
   UpdateScoring();

   // if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
}

G4bool G4tbbWorkerRunManager::DoOneEvent( G4int i_event, G4String msg)
{
   // This method is called by the tbb::task. 
   // It 
   //   i) performs the simulation of a single event
   //  ii) obtains the status
   // iii) stacks the event (!?)

   // Set the UI command text 
   msgText= msg;

   // Do the basic work
   ProcessOneEvent( i_event ); 

   StackPreviousEvent(currentEvent);
   currentEvent = 0;
   //G4cout<<"runAborted is:"<<runAborted<<G4endl;
   return runAborted;
}

void G4tbbWorkerRunManager::DoEventLoop(G4int ,      // n_event, 
                                        const char*, // macroFile , 
                                        G4int )      // n_select)
{
   G4cerr << " ERROR> The event Loop is not handled by the (worker) "
          << " Run Manager - it is handled by TBB" << G4endl;
}



G4tbbWorkerRunManager::G4tbbWorkerRunManagerInstancesType& 
  G4tbbWorkerRunManager::GetInstancesList()
{
  return instancesList; 
}

unsigned int G4tbbWorkerRunManager::NumberOfWorkers()
{
   return instancesList.unsafe_size();
}

void G4tbbWorkerRunManager::DestroyWorkersAndCleanup()
{
   // instancesList.ClearAndDestroy(); 
}

