#include "ParN02Job.hh"
#include <G4tbbRunManager.hh>
#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include <G4VisExecutive.hh>
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "FTFP_BERT.hh"

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

//Called once by main thread
void ParN02Job::UserActions(G4tbbRunManager* /*rm*/) 
{
  G4cout<<"ParN02Job create detector, start of UserActions"<<G4endl;
  detector = new ExN02DetectorConstruction();
  G4cout<<"ParN02Job detector created, end of UserActions"<<G4endl;
}

//Called by each working thread
void ParN02Job::InitSetup(G4tbbRunManager* )  // rm
{
  G4cout<<"ParN02Job InitSetup : Worker Consutrction of SD and field"<<G4endl;
  assert( detector );
  // detector->SlaveExN02DetectorConstruction();
  detector->ConstructSDandField();
  G4cout<<"ParN02Job InitSetup : Worker Construction - done"<<G4endl;
}

//This is common between threads, basically a copy of the original main part...
void ParN02Job::JobPrepare(G4tbbRunManager* rm )
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
  
// #ifdef G4VIS_USE
#if 0
  // Visualization, if you choose to have it!
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  G4VUserPrimaryGeneratorAction* gen_action = 
     new ExN02PrimaryGeneratorAction(detector);
  rm->SetUserAction(gen_action);

  G4UserRunAction* run_action = new ExN02RunAction;
  rm->SetUserAction(run_action);  

  G4UserEventAction* event_action = new ExN02EventAction;
  rm->SetUserAction(event_action);

  G4UserSteppingAction* stepping_action = new ExN02SteppingAction;
  rm->SetUserAction(stepping_action);
  G4cout<<"ParN02Job JobPrepare : done"<<G4endl;

}
