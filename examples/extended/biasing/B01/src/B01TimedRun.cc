#include "B01TimedRun.hh"
#include "G4RunManager.hh"
#include "B01PhysicsList.hh"
#include "B01PrimaryGeneratorAction.hh"
#include "B01DetectorConstruction.hh"
#include "B01EventAction.hh"
#include "ExN03RunAction.hh"
#include "G4CellScorer.hh"


//#include "G4DosimPhysics.hh"
//#include "B01SteppingAction.hh"

B01TimedRun::B01TimedRun(G4int time)
  :
  fRunManager(new G4RunManager),
  fTime(time)
{
  if (!fRunManager) {
    G4std::G4Exception("B01TimedRun::B01TimedRun: new failed to create G4RunManager!");
  }
}
B01TimedRun::~B01TimedRun(){}

void B01TimedRun::BeamOn(G4int n){
  if (!fRunManager) {
    G4std::G4Exception("B01TimedRun::BeamOn: no G4RunManager");
  }
  fRunManager->BeamOn(n);
}


void B01TimedRun::SetSpecialG4CellScorer(const G4CellScorer *cellscorer){
  fRunManager->SetUserAction(new B01EventAction(cellscorer, fTime));
}

void B01TimedRun::SetDetector(G4VPhysicalVolume &worldvol){

  
  fRunManager->
    SetUserInitialization(new B01DetectorConstruction(worldvol));
  
  fRunManager->SetUserInitialization(new B01PhysicsList);


  //  fRunManager->SetUserInitialization(new PhysicsList);
  //  fRunManager->SetUserAction(new B01SteppingAction);



  fRunManager->SetUserAction(new B01PrimaryGeneratorAction);  

}
void B01TimedRun::Initialize(){
  
  fRunManager->Initialize();
}


