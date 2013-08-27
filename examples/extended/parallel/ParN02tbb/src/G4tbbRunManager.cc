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
#include "G4Timer.hh"
#include "G4tbbRunManager.hh"
#include "G4tbbTask.hh"
#include "G4UImanager.hh"
#include "G4VtbbJob.hh"
#include "tbbhelper.hh"

#include "G4tbbWorkerRunManager.hh"

//TBB includes
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>

//Global instance of seeds queue
G4tbbRunManager::G4tbbSeedsQueueType 
  G4tbbRunManager::seedsQueue = G4tbbRunManager::G4tbbSeedsQueueType();

G4tbbRunManager* G4tbbRunManager::fMasterTbbRM = 0; 

G4tbbRunManager::G4tbbRunManager() : G4RunManager() ,
  seedsSequenceLength(1),
  fFirstBeamOnCall(true),
  fCreatedTaskList(false),
  fNumWorkers(1)
{
   if( fMasterTbbRM ) 
        G4Exception("G4RunManager", "Run034",
                    FatalException, "Must only define one *master* TBB Run Manager." ); 
   else
      fMasterTbbRM= this; 
}

G4tbbRunManager::~G4tbbRunManager()
{
  // Destroy the Worker run managers created to carry out the tasks
  G4tbbWorkerRunManager::DestroyWorkersAndCleanup(); 

  //  These should call the cleanup methods on each worker (thread) used
  //    the methods will be similar to 
  //      tbbSlaveDestroyGeometryAndPhysicsVector();
  
  if( fCreatedTaskList ) delete fpTaskList;
}

void G4tbbRunManager::InitializeAndCreateTaskList( G4int noWorkers )
{
  G4cout<<"Initializing task_manager with "<< noWorkers<<" workers."<<G4endl;
  tbb::task_scheduler_init init( noWorkers );
  
  fNumWorkers= noWorkers;
  
  G4cout<<"Creating TBB tasklist and setting it to the runManager"<<G4endl;
  tbb::task_list *tlist= new tbb::task_list;
  SetTaskList( tlist );
  assert(tlist->empty()==true);
  
  fCreatedTaskList=true;
}


void G4tbbRunManager::BeamOn(G4int n_event, const char* macroFile,G4int num_select)
{
  
  G4cout<<"BeamOn called. ";
  G4cout<<" Thread ID is:"<< gettid()
  <<"; runManager pointer:"<<this<<G4endl;
  
  if(!fpTaskList) InitializeAndCreateTaskList( fNumWorkers );
  
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
  {
    numberOfEventToBeProcessed = n_event;
    ConstructScoringWorlds();
    RunInitialization();
    if(n_event>0)
    {
      DoEventLoop(n_event,macroFile,num_select);
      assert( fpTaskList->empty()==false );
    }
  }
  
  SpawnTasksAndWait(fpTaskList);
  
  G4cout<<"tbbRunManager::BeamOn> Work done"<<G4endl;
  
  G4cout<<"tbbRunManager::BeamOn> Calling RunTermination()."<<G4endl;
  RunTermination();

  // Code copied from standard RunManager method TerminateEventLoop()
  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    {
      G4cout << "  Run Aborted " << G4endl;
      // "after " << i_event + 1  << " events processed." << G4endl;
    }
    else
    {
      G4cout << "  Number of events processed : " << n_event << G4endl;
    }
    G4cout << "  "  << *timer << G4endl;
  }
  
  fFirstBeamOnCall= false;
}


void G4tbbRunManager::DoEventLoop( G4int n_event,
                                  const char* macroFile,
                                  G4int n_sel) {
  
  //Implementation of G4 event loop
  //Problem : this must be done only once (per run).
  //  This means that the BeamOn command must be applied only once,
  //   ie by the "master thread" (in the sense of G4MT)
  //  if ( isSlave != 0 ) return;
  
  n_select = n_sel;
  //Copy from standard RunManager
  if(verboseLevel>0)
  { timer->Start(); }
  
  if(macroFile!=0)
  {
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }
  
  // "Event loop", in this case it
  //   Creates a list of tasks: one event= one task
  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ )
  {
    G4tbbTask& task = *new(tbb::task::allocate_root())
    G4tbbTask(job,i_event,/*n_select,msg,*/seedsSequenceLength);
    fpTaskList->push_back( task );
  }
  G4cout<<"G4tbbRunManager::DoEventLoop> Created " << n_event << " tasks."<<G4endl;
  
}

void G4tbbRunManager::SpawnTasksAndWait( tbb::task_list* tlist )
{
  assert(!fpTaskList);
  
  try {
    G4cout<<"Now spawn work and waiting"<<G4endl;
    tbb::task::spawn_root_and_wait( *tlist );
  } catch(std::exception& e) {
    G4cerr<<"Error occurred. Error test is:\""<<e.what()<<"\""<<G4endl;
  }
}

void G4tbbRunManager::AddSeed( const long& seed ) {
  seedsQueue.push( seed );
}

G4bool G4tbbRunManager::GetNextSeed( long& seed ) {
  return seedsQueue.try_pop( seed );
}

unsigned int G4tbbRunManager::GetNumberOfWorkers()
{
  return G4tbbWorkerRunManager::NumberOfWorkers();
}

/******
G4bool G4tbbRunManager::DoOneEvent( G4int i_event)
       //, G4int n_select, const G4String& msg)
{
  // This method is called by the tbb::task. 
  // It perfomrs the simulation of a single event
  currentEvent = GenerateEvent(i_event);
  eventManager->ProcessOneEvent(currentEvent);
  AnalyzeEvent(currentEvent);
  UpdateScoring();
  if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msg);
  StackPreviousEvent(currentEvent);
  currentEvent = 0;
  //G4cout<<"runAborted is:"<<runAborted<<G4endl;
  return runAborted;
}
 ******/ 
