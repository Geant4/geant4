#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include <unistd.h>
#define _GNU_SOURCE
#include <getopt.h>
#include "g4std/set"
#include "G4UImanager.hh"


#include <iomanip>

#include "B03DetectorConstruction.hh"
#include "B03PhysicsList.hh"
#include "B03PrimaryGeneratorAction.hh"

// Files specific for biasing
#include "G4MassImportanceManager.hh"


int main(int argc, char **argv) {
  

  int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the detector      ---------------------------
  B03DetectorConstruction *detector = new B03DetectorConstruction();
  runManager->
    SetUserInitialization(detector);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B03PhysicsList);
  runManager->SetUserAction(new B03PrimaryGeneratorAction);
  runManager->Initialize();

  // the IStore is filled during detector construction
  G4VIStore &aIstore = *detector->GetIStore();
  // create the importance manager for biasing in the tracking world
  G4MassImportanceManager mim(aIstore, "neutron");
  mim.Initialize();


  G4UImanager* UI;

  UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/control/execute init.mac");   


  runManager->BeamOn(numberOfEvent);

  return 0;
}

