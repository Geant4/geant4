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
/*
Author: Susanna Guatelli
*/
//
//    **********************************
//    *                                *
//    *     BrachyPhysicsList.cc       *
//    *                                *
//    **********************************
//
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "BrachyPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyPhysicsListMessenger.hh"

BrachyPhysicsList::BrachyPhysicsList():  G4VModularPhysicsList()
{
SetVerboseLevel(1); 
 
// EM physics: default
fEmPhysicsList = new G4EmLivermorePhysics();
fEmName="emlivermore";

// Add Decay
fDecPhysicsList = new G4DecayPhysics();
fRadDecayPhysicsList = new G4RadioactiveDecayPhysics();
fMessenger = new BrachyPhysicsListMessenger(this);
}

BrachyPhysicsList::~BrachyPhysicsList()
{  
delete fMessenger;
delete fDecPhysicsList;
delete fRadDecayPhysicsList;
delete fEmPhysicsList;
}

void BrachyPhysicsList::ConstructParticle()
{
fDecPhysicsList -> ConstructParticle();
}

void BrachyPhysicsList::ConstructProcess()
{
AddTransportation();
fEmPhysicsList -> ConstructProcess();

// decay physics list
fDecPhysicsList -> ConstructProcess();
fRadDecayPhysicsList -> ConstructProcess();
}


void BrachyPhysicsList::AddPhysicsList(const G4String& name)
{
  
  if (name == fEmName) return;

  if (name == "emstandard_opt0"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  } else if (name == "emstandard_opt1"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt2"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt3"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt4"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  } else if (name == "empenelope"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  } else if (name == "emlivermore"){
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;  
  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
  G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is activated"
           << G4endl;
}

void BrachyPhysicsList::SetCuts()
{
// Definition of  threshold of production 
// of secondary particles
// This is defined in range.
defaultCutValue = 0.05 * mm;
SetCutValue(defaultCutValue, "gamma");
SetCutValue(defaultCutValue, "e-");
SetCutValue(defaultCutValue, "e+");
  
// By default the low energy limit to produce 
// secondary particles is 990 eV.
// This value is correct when using the EM Standard Physics.
// When using the Low Energy Livermore this value can be 
// changed to 250 eV corresponding to the limit
// of validity of the physics models.
// Comment out following three lines if the 
// Standard electromagnetic Package is adopted.
G4double lowLimit = 250. * eV;
G4double highLimit = 100. * GeV;

G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit,
                                                                highLimit);

// Print the cuts 
if (verboseLevel>0) DumpCutValuesTable();
}
