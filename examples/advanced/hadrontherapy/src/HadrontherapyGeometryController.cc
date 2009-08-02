#include "HadrontherapyGeometryController.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "IAEADetectorConstruction.hh"
#include "G4RunManager.hh"

HadrontherapyGeometryController::HadrontherapyGeometryController()
{

}

HadrontherapyGeometryController::~HadrontherapyGeometryController()
{

}

void HadrontherapyGeometryController::SetGeometry(G4String name)
{
  G4cout <<"Activating geometry " << name << G4endl;
  if(name == "IAEA") {
    registerGeometry(new IAEADetectorConstruction());
    G4cout <<"IAEA geometry activated" << G4endl;
  } else if(name == "default") {
    registerGeometry(new HadrontherapyDetectorConstruction());
  } else {
    G4cout <<"Unknown geometry: " << name << ". Geometry not changed." << G4endl;
  }
}

void HadrontherapyGeometryController::registerGeometry(G4VUserDetectorConstruction *detector)
{
  G4RunManager *runManager = G4RunManager::GetRunManager();
  runManager->SetUserInitialization(detector);
  runManager->GeometryHasBeenModified();
}

