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
// Previous authors: G. Guerrieri and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
// 
//
//    **********************************
//    *                                *
//    *     G4HumanPhantomPhysicsList.cc       *
//    *                                *
//    **********************************
//
//
#include "G4HumanPhantomPhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"              
#include "G4VModularPhysicsList.hh" 
#include "G4EmStandardPhysics_option4.hh" 
#include "G4DecayPhysics.hh"

G4HumanPhantomPhysicsList::G4HumanPhantomPhysicsList():  G4VModularPhysicsList()
{
  SetVerboseLevel(1);
  emPhysicsList = new G4EmStandardPhysics_option4(1);
  decPhysicsList = new G4DecayPhysics();
// Alternatively you can substitute this physics list
// with the LowEnergy Livermore or LowEnergy Penelope: 
// emPhysicsList = new G4EmLivermorePhysics();
// Low Energy based on Livermore Evaluated Data Libraries
//
// Penelope physics
//emPhysicsList = new G4EmPenelopePhysics();
}

G4HumanPhantomPhysicsList::~G4HumanPhantomPhysicsList()
{
delete decPhysicsList;
delete emPhysicsList;
}

void G4HumanPhantomPhysicsList::ConstructParticle()
{
 decPhysicsList -> ConstructParticle(); 
}

void G4HumanPhantomPhysicsList::ConstructProcess()
{
  AddTransportation();
  emPhysicsList -> ConstructProcess();
  decPhysicsList -> ConstructProcess();
}

void G4HumanPhantomPhysicsList::SetCuts()
{
  // The production threshold is fixed to 0.1 mm for all the particles
  // Secondary particles with a range bigger than 0.1 mm 
  // are generated; otherwise their energy is considered deposited locally

  defaultCutValue = 1. * mm;

  const G4double cutForGamma = defaultCutValue;
  const G4double cutForElectron = defaultCutValue;
  const G4double cutForPositron = defaultCutValue;

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  // Set the secondary production cut lower than 990. eV
  // Very important for high precision of lowenergy processes at low energies
 
  G4double lowLimit = 250. * eV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
  
  if (verboseLevel>0) DumpCutValuesTable();
}
