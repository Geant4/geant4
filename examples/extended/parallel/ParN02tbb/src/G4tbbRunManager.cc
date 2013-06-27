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

//Global instance of seeds queue
G4tbbRunManager::G4tbbSeedsQueueType 
  G4tbbRunManager::seedsQueue = G4tbbRunManager::G4tbbSeedsQueueType();

G4tbbRunManager::G4tbbRunManagerInstancesType 
 G4tbbRunManager::instancesList=G4tbbRunManager::G4tbbRunManagerInstancesType();

G4tbbRunManager::G4tbbRunManagerInstancesType& 
  G4tbbRunManager::GetInstancesList()
{
  return instancesList; 
}

// Static TLS member if needs to perform intialization of "slave thread" 
//  in the sense of G4MT i.e. content of my_slave_thread
//tbb::enumerable_thread_specific< G4bool > G4tbbRunManager::isSlaveInitialized 
//  = false;


G4tbbRunManager::G4tbbRunManager() : G4RunManager() ,
  seedsSequenceLength(1)
{
  instancesList.push(this);
  //if ( isSlave !=0 && isSlaveInitialized == false ){
  //  InitSlave();
  //}
}

G4tbbRunManager::G4tbbRunManager( int isSlaveFlag ) : G4RunManager(isSlaveFlag),
  seedsSequenceLength(1)
{
  instancesList.push(this);
  //if ( isSlave !=0 && isSlaveInitialized == false ){
  //  InitSlave();
  //}
}

G4tbbRunManager::~G4tbbRunManager()
{
  //if ( isSlave != 0 ) 
  //  tbbSlaveDestroyGeometryAndPhysicsVector(); 
  //    ^^^^ How to call this on each thread???
}

//void G4tbbRunManager::InitSlave(){
  //Here do what is done in my_slave_thread, namely:
  //tbbSlaveBuildGeometryAndPhysicsVector();
  //G4_main() <<<---- Ah Ah!
//}
void G4tbbRunManager::DoEventLoop( G4int n_event,
                                   const char* macroFile,
                                   G4int n_sel) {

  //Implementation of G4 event loop
  //Problem : this should be done only once. This means that the BeamOn 
  // command should be applied only once, 
  //   ie by the "master thread" (in the sense of G4MT)
  if ( isSlave != 0 ) return;
  
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

  //Event loop, in this case it's a list of tasks: one event= one task
  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ )
  {
    G4tbbTask& task = *new(tbb::task::allocate_root()) 
               G4tbbTask(job,i_event,/*n_select,msg,*/seedsSequenceLength);
    tasklist->push_back( task );
  }
  
  // Code copied from standard RunManager method TerminateEventLoop()
  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { 
       G4cout << "  Run Aborted after " << i_event + 1 
              << " events processed." << G4endl; 
    }
    else
    { 
       G4cout << "  Number of events processed : " << n_event << G4endl; 
    }
    G4cout << "  "  << *timer << G4endl;
  }
}

void G4tbbRunManager::AddSeed( const long& seed ) {
  seedsQueue.push( seed );
}

G4bool G4tbbRunManager::GetNextSeed( long& seed ) {
  return seedsQueue.try_pop( seed );
}

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
