#include "B03App.hh"
#include "g4std/iostream"

#include "CLHEP/Random/Random.h"

#include "G4RunManager.hh"

#include "G4VPhysicalVolume.hh"
#include "B03DetectorConstruction.hh"
#include "B03PhysicsList.hh"
#include "B03PrimaryGeneratorAction.hh"


B03AppBase *B03AppBase::fB03AppBase = 0;
B03AppBase &B03AppBase::GetB03AppBase() {
  if (!fB03AppBase) {
    fB03AppBase = new B03AppBase;
  }
  return *fB03AppBase;
}

B03AppBase::B03AppBase(){
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  frunMgr = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  fDetector = new B03DetectorConstruction;
  frunMgr->SetUserInitialization(fDetector);
  //  ---------------------------------------------------
  fPhysics = new B03PhysicsList;
  frunMgr->SetUserInitialization(fPhysics);
  fPrimary = new B03PrimaryGeneratorAction;
  frunMgr->SetUserAction(fPrimary);
  frunMgr->Initialize();

  /*
  // create the detector for the importances and scoring
  fImportancedetector
    = new B03ImportanceDetectorConstruction();
  // the IStore is filled during detector construction

  fIStore = fImportancedetector->GetIStore();
  
    // create importance sampler to importance sample neutrons
    // in the importance geometry
  fCS_store = 
    &fImportancedetector->GetCellScorerStore();
  fCS_scorer = new G4CellStoreScorer(*fCS_store);
  fPmgr = new
    G4ParallelImportanceScoreSampler(*fIStore, *fCS_scorer, "neutron");
  fPmgr->Initialize();
  */
}

B03AppBase::~B03AppBase(){
  delete fPrimary;
  delete fPhysics;
  delete fDetector;
  delete frunMgr;
}


