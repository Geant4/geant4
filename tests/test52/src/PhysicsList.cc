
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "Particles.hh"
#include "EMElectronEEDL.hh"
#include "EMElectronPenelope.hh"
#include "EMElectronStandard.hh"
#include "EMPositronPenelope.hh"
#include "EMPositronStandard.hh"
#include "EMPhotonEPDL.hh"
#include "EMPhotonPenelope.hh"
#include "EMPhotonStandard.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

G4VPhysicsConstructor* PhysicsList::emElectron = 0;
G4VPhysicsConstructor* PhysicsList::emPositron = 0;
G4VPhysicsConstructor* PhysicsList::emPhoton = 0;


PhysicsList::PhysicsList() :
    G4VModularPhysicsList() {

  SetVerboseLevel(1);
  defaultCutValue = 0.001 * mm;   

  messenger = new PhysicsListMessenger(this);

  RegisterPhysics(new Particles());
}


PhysicsList::~PhysicsList() {

  delete messenger;
}


void PhysicsList::RegisterPhysConstructor(const G4String& constrName) {

  // ****************
  // *** Electron ***
  // ****************
  
  if(emElectron == 0) {
     if(constrName == "Electron-EEDL") 
        emElectron = new EMElectronEEDL();
     if(constrName == "Electron-Penelope") 
        emElectron = new EMElectronPenelope();
     if(constrName == "Electron-Standard") 
        emElectron = new EMElectronStandard();

     if(emElectron) RegisterPhysics(emElectron);
  }

  // ****************
  // *** Positron ***
  // ****************

  if(emPositron == 0) {
     if(constrName == "Positron-Penelope") 
        emPositron = new EMPositronPenelope();
     if(constrName == "Positron-Standard") 
        emPositron = new EMPositronStandard();

     if(emPositron) RegisterPhysics(emPositron);
  }

  // **************
  // *** Photon ***
  // **************

  if(emPhoton == 0) {
     if(constrName == "Photon-EPDL") 
        emPhoton = new EMPhotonEPDL();
     if(constrName == "Photon-Penelope") 
        emPhoton = new EMPhotonPenelope();
     if(constrName == "Photon-Standard") 
        emPhoton = new EMPhotonStandard();

     if(emPhoton) RegisterPhysics(emPhoton);
  }
}


void PhysicsList::SetProdThreshold(double cut) {

  if(cut > 0.0 * mm) {
     defaultCutValue = cut;
  }
}


void PhysicsList::SetCuts() {

  G4double lowerProdLimit = 250.0 * eV;
  G4double upperProdLimit = 1.0 * GeV;

  G4ProductionCutsTable::GetProductionCutsTable()
      -> SetEnergyRange(lowerProdLimit, upperProdLimit);

  SetCutsWithDefault();

  G4Region* targetRegion = 
      G4RegionStore::GetInstance() -> GetRegion("Target");

  G4ProductionCuts* prodCutsTarget = new G4ProductionCuts;
  prodCutsTarget -> SetProductionCut(defaultCutValue);

  targetRegion -> SetProductionCuts(prodCutsTarget);

  DumpCutValuesTable();
}
