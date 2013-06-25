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
// $Id$
// GEANT4 tag $Name: geant4-09-00 $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#include "ExN02DetectorConstruction.hh"
//#include "ExN02PhysicsList.hh"
//#include "ExN02PrimaryGeneratorAction.hh"
//#include "ExN02RunAction.hh"
//#include "ExN02EventAction.hh"
//#include "ExN02SteppingAction.hh"
#include "ParN02Job.hh"

#include "ExN02SteppingVerbose.hh"

#include "G4tbbRunManager.hh"

//#include "G4RunManager.hh"
//#include "G4UImanager.hh"

//#ifdef G4VIS_USE
//#include "G4VisExecutive.hh"
//#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//01.25.2009 Xin Dong: This example came from the original sequential
//program FullCMS. The original program is changed here to support parallel
//computing with multiple threads. All events are assigned to each worker
//thread in a round robin fashion. All threads share most detector data 
//including physics table and physics vector for some physics processes.
//The master process initializes the data in a regular way. However, worker
//threads initialize thread private data only.
//#include "G4MTParTopC.icc"
void my_slave_thread(void*) {} 
//01.25.2009 Xin Dong: Threads share this object.
//ExN02DetectorConstruction* detector = 0;

#include <CLHEP/Random/RanluxEngine.h>

//TBB includes
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>
// #include "/Users/japost/Software/TBB/4.1/include/tbb/task_scheduler_init.h"
// #include "/Users/japost/Software/TBB/4.1/include/tbb/task.h"

#include "G4MTGetTid.hh"

int main(int argc,char** argv)
{
  assert( G4RunManager::GetRunManager() == NULL );
  //Number of events
  int nevents = 10;
  //Macro file
  G4String aMacroFileName = "run.mac";
  //Number of threads
  int nthreads = 2;
  for ( int idx = 1 ; idx<argc ; ++idx ) {
    if ( strcmp(argv[idx],"-n")==0 ) nevents = atoi( argv[idx+1] );
    if ( strcmp(argv[idx],"-t")==0 ) nthreads =atoi( argv[idx+1] );
    if ( strcmp(argv[idx],"-m")==0 ) aMacroFileName = argv[idx+1];
  }
  //Random engine
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  //  CLHEP::HepRandom::setTheSeed( seed ); 
  CLHEP::HepRandom::setTheSeed( 1220515164 );
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
         << " ===================================================== " << G4endl 
         << G4endl; 
  G4cout<<"Random Engine: "<<CLHEP::HepRandom::getTheEngine()<<G4endl;
  //A G4 stuff...
  G4VSteppingVerbose* verbosity = new ExN02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  //Create a G4VtbbJob object: it contains all job stuff
  G4VtbbJob* aJob = new ParN02Job( aMacroFileName );

  //Create the G4tbbRunManager object, and set the associated job
  G4cout<<"Creating G4tbbRunManager"<<G4endl;
  G4tbbRunManager* runManager = new G4tbbRunManager();
  runManager->SetJob( aJob );
  //Set the list of seeds. Depending on the engine type, set the correct 
  // number of seeds per event. In this case 1 seed / evt.
  G4cout<<"Adding random seeds to the runManager"<<G4endl;
  unsigned int aseed;
  for ( int i = 0 ; i < nevents ; ++i ) {
    aseed = static_cast<unsigned int>(*CLHEP::HepRandom::getTheEngine());
    G4cout<<"Seed number "<<i<<" : "<<aseed<<G4endl;
    runManager->AddSeed( aseed );
  }

  //Now create TBB stuff...
  //List of tasks
  G4cout<<"Initializing task_manager with "<<nthreads<<" threads."<<G4endl;
  tbb::task_scheduler_init init( nthreads );

  G4cout<<"Creating TBB tasklist and setting it to the runManager"<<G4endl;
  tbb::task_list tasklist;
  runManager->SetTaskList( &tasklist );
  assert(tasklist.empty()==true);

  //Now start....
  G4cout<<"Calling G4VtbbJob::InitRun"<<G4endl;
  aJob->InitRun( runManager );
  G4cout<<"Now calling beamOn"<<G4endl;
  runManager->BeamOn(nevents);
  assert( tasklist.empty()==false );
  
  G4cout<<"BeamOn called"<<G4endl;
  G4cout<<"My Thread ID is:"<<gettid()
       <<"; runManager pointer:"<<runManager<<G4endl;
  try {
    G4cout<<"Now spawn work and waiting"<<G4endl;
    tbb::task::spawn_root_and_wait( tasklist );
  } catch(std::exception& e) {
    G4cerr<<"Error occurred. Error test is:\""<<e.what()<<"\""<<G4endl;
  }
  G4cout<<"Work done"<<G4endl;

  G4cout<<G4endl<<"End of job, deleting stuff"<<G4endl;
  delete aJob;
  //delete runManager;
  //Delete all instances of runManger
  G4tbbRunManager::G4tbbRunManagerInstancesType queue = 
     G4tbbRunManager::GetInstancesList();
  G4cout<<"Number of G4tbbRunManager instances: "<<queue.unsafe_size()<<G4endl;
  G4tbbRunManager* ist = 0;
  //while ( queue.try_pop( ist ) ) { delete ist;}
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

