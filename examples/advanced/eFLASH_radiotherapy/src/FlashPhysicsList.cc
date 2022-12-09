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
//
/// \file FlashPhysicsList.cc
/// \brief Implementation of the FlashPhysicsList class

#include "FlashPhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4ProductionCuts.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"

FlashPhysicsList::FlashPhysicsList() : G4VModularPhysicsList() {
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

 

  // EM physics
  //RegisterPhysics(new G4EmPenelopePhysics());
  RegisterPhysics(new G4EmStandardPhysics_option4());
}

FlashPhysicsList::~FlashPhysicsList() {}

void FlashPhysicsList::SetCuts() {//set cuts based on region name defined in detector construction

  SetCutsWithDefault();
  G4Region *region;
  G4String regName;
  G4ProductionCuts *cuts;



  regName = "Phantom_reg";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("e+"));
  region->SetProductionCuts(cuts);
  
 
}

