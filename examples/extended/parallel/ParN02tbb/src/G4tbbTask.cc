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
#include "G4tbbTask.hh"
#include "G4tbbRunManager.hh"
//#define AUTOLOCKDEBUG
#include "G4AutoLock.hh"
#include "G4VtbbJob.hh"
#include "tbbhelper.hh"

G4tbbTask::G4tbbTask(G4VtbbJob* j , G4int ev, G4int l) :
event(ev), job(j)

{ 
  //Note reprorducibility is guaranteed if G4tbbTask are created
  //sequentially, in this case the event number event will "pop"
  //always the same seeds.
  //The seeds is a zero terminated array.
  seeds = new long[l + 1];
  for ( int i = 0 ; i<l ; ++i )
    G4tbbRunManager::GetNextSeed( *(seeds + i) );
  seeds[l] = 0;

  //  Extension: add G4int sel, const G4String& m as arguments,
  //  and create before { :
  //       select(sel),  msg(m)
}

G4tbbTask::~G4tbbTask() {
  delete seeds;
}


G4Mutex  runmanagermutex = G4MUTEX_INITIALIZER;
#include <tbb/enumerable_thread_specific.h>

typedef tbb::enumerable_thread_specific<G4bool> InitType;
InitType isInitialized;
InitType canGo;
#include <tbb/compat/thread>

tbb::task* G4tbbTask::execute() {
  G4RunManager* runmanager = G4RunManager::GetRunManager();
  G4tbbRunManager* rm = 0;
  TBBMSG("Starting execute. Event number: "<<event
         <<" RunManager pointer:"<<runmanager);
  if ( runmanager == NULL ) {
    //wasting some time (or the fast thread will use all events)
    //std::this_thread::sleep_for( tbb::tick_count::interval_t(0.05) );
    //return NULL; //skip this event
    TBBMSG("No valid runmanager found, creating one.");
    {    
       G4AutoLock lock( &runmanagermutex );

       //If I am here it means that this is a new thread.
       //  I need to prepare the G4 state for the thread
       //    and then do the needed work

       //Step 1- Copy geometry and physics tables
       tbbSlaveBuildGeometryAndPhysicsVector();

       //Step 2- Create a new (worker) run manager
       G4int isSlave = 1;
       G4iosInitialization();//This is needed or next line will crash!
       rm = new G4tbbRunManager( isSlave );
       TBBMSG("End of creation of slave runmanager");
       // lock.explicit_unlock();
    }
    //Step 3- Perform initialization
    TBBMSG("Calling threadsafeinitsetup");
    job->ThreadSafeInitSetup( rm );
    //return NULL;////AAAAA
  } else { 
    rm = static_cast<G4tbbRunManager*>( runmanager );
  }
  assert(rm!=NULL);
  //First thing to do: set the seeds, for next event
  //G4cout<<"In Execute the random engine pointer is: "
  //      <<CLHEP::HepRandom::getTheEngine()<<G4endl;
  CLHEP::HepRandom::setTheSeeds( seeds , -1 );
  //Now do what is done in BeamOn to initialize the run...
  InitType::reference myInit = isInitialized.local();
  InitType::reference myCanGo= canGo.local();
  if ( myInit == false ) {
    myCanGo = true;
    TBBMSG("This only once");
    //Need to do stuff...
    G4bool cond = rm->ConfirmBeamOnCondition();
    if ( cond ) {
      rm->ConstructScoringWorlds();
      rm->RunInitialization();//From now on ConfirmBeamOnCondition == false
    }
    else {
      myCanGo = false;
    }
    myInit = true;
  }
  G4bool abortrun = false;
  if ( myCanGo ) // Check that initialization has been done correctly
    abortrun = rm->DoOneEvent( event);
  
  if ( abortrun == true ) {
    G4cout<<"ERROR:: ABORTRUN REQUEST!!!!"<<G4endl;
  }
  //How to do RunTermination at the end of the run????
  return static_cast<tbb::task*>(NULL);
 }
