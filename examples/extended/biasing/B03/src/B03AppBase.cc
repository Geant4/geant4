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
// $Id: B03AppBase.cc,v 1.3 2006/06/29 16:35:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// --------------------------------------------------------------------
#include "B03AppBase.hh"
#include <iostream>

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

  CLHEP::HepRandom::setTheSeed(myseed);

  frunMgr = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  fDetector = new B03DetectorConstruction;
  frunMgr->SetUserInitialization(fDetector);
  //  ---------------------------------------------------
  fPhysics = new B03PhysicsList;
  frunMgr->SetUserInitialization(fPhysics);
  fPrimary = new B03PrimaryGeneratorAction;
  frunMgr->SetUserAction(fPrimary);
  //  frunMgr->Initialize();

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


