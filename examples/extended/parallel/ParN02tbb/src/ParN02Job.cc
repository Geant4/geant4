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
#include "ParN02Job.hh"
#include "G4WorkerRunManager.hh"
#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "FTFP_BERT.hh"

#include "G4VisExecutive.hh"

ParN02Job::ParN02Job(const G4String& mf) :
G4VtbbJob(mf),
detector(0)
{
  G4cout<<"ParN02Job constructor. Executing macro: "<<mf<<G4endl;
}

ParN02Job::~ParN02Job() {
  G4cout<<"ParN02Job destructor"<<G4endl;
}

//Called once by main thread
void ParN02Job::CreateDetector(G4tbbRunManager* /*rm*/)
{
  G4cout<<"ParN02Job Create detector, start."<<G4endl;
  detector = new ExN02DetectorConstruction();
  G4cout<<"ParN02Job detector created,  end."<<G4endl;
}

//Called by worker thread - which are not the master 
void ParN02Job::InitWorkerSetup(G4tbbWorkerRunManager* )  // rm
{
  G4cout<<"ParN02Job InitSetup : Worker Consutrction of SD and field"<<G4endl;
  assert( detector );
  detector->ConstructSDandField();
  G4cout<<"ParN02Job InitSetup : Worker Construction - done"<<G4endl;
}

// This is common between threads, basically a copy of the 
//   'old' main function - with only the User Action creation 
//   separated.
//
void ParN02Job::JobPrepare(G4RunManager* rm )
{
  //G4cout<<"PIPPO Random:"<<G4Random::getTheEngine()<<G4endl;
  //These two guarantee to use the same random generator for all threads
  //CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  //G4Random::setTheEngine( &defaultEngine ); 
  
  G4cout<<"ParN02Job JobPrepare : start"<<G4endl;

  rm->SetUserInitialization(detector);
  G4VUserPhysicsList* physics = new FTFP_BERT;//ExN02PhysicsList;
  rm->SetUserInitialization(physics);

  // Altenative: Obtain (or create) a "Physics Workspace"
  //
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  G4cout<<"ParN02Job JobPrepare : done"<<G4endl;
}

//Called once by each worker thread - this can include the master
void ParN02Job::UserActions(G4RunManager* rm)
{
  G4cout<<"ParN02Job: start of UserActions"<<G4endl;
  G4VUserPrimaryGeneratorAction* gen_action = 
     new ExN02PrimaryGeneratorAction(detector);
  rm->SetUserAction(gen_action);

  G4UserRunAction* run_action = new ExN02RunAction;
  rm->SetUserAction(run_action);  

  G4UserEventAction* event_action = new ExN02EventAction;
  rm->SetUserAction(event_action);

  G4UserSteppingAction* stepping_action = new ExN02SteppingAction;
  rm->SetUserAction(stepping_action);

  G4cout<<"ParN02Job: end of UserActions"<<G4endl;
}
