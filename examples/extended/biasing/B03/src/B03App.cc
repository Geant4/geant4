//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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


