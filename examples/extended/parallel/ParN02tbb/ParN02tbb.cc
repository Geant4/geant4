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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void my_slave_thread(void*) {} 

#include "Randomize.hh"

//TBB includes
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>

#include "G4MTGetTid.hh"

int main(int argc,char** argv)
{
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
  G4Random::setTheEngine( &defaultEngine ); 

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
  unsigned long aseed;

  for ( int i = 0 ; i < 2*nevents ; ++i ) {
    aseed = (unsigned long) (100000000L * G4Random::getTheGenerator()->flat());
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

  G4cout<< "Number of Tbb Worker instances used: : "
        <<  runManager->GetNumberOfWorkers() <<G4endl;
  // The TBB runManager 'controls' all its worker Run Managers - it will delete them
  delete runManager;

  // G4tbbRunManager* ist = 0;
  // while ( queue.try_pop( ist ) ) { delete ist;}
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

