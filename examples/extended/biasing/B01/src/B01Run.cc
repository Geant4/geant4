#include "B01Run.hh"
#include "G4RunManager.hh"
#include "B01PhysicsList.hh"
#include "B01PrimaryGeneratorAction.hh"
#include "B01DetectorConstruction.hh"


B01Run::B01Run()
  :
  fRunManager(new G4RunManager)
{}
B01Run::~B01Run(){}
void B01Run::SetDetector(G4VPhysicalVolume &worldvol){

  
  fRunManager->
    SetUserInitialization(new B01DetectorConstruction(worldvol));
  
  fRunManager->SetUserInitialization(new B01PhysicsList);
  fRunManager->SetUserAction(new B01PrimaryGeneratorAction);

}
void B01Run::Initialize(){
  fRunManager->Initialize();
}

void B01Run::BeamOn(G4int nevents) const {
  fRunManager->BeamOn(nevents);
}

