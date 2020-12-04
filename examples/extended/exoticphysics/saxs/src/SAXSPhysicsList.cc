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
/// \file SAXSPhysicsList.cc
/// \brief Implementation of the SAXSPhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SAXSPhysicsList.hh"
#include "SAXSPhysicsListMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4ProcessManager.hh"
#include "G4Threading.hh"

//Particles
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//EM Physics Lists
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmPenelopePhysicsMI.hh"

//EM options
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

//Hadronic and Extra Physics Lists
#include "G4EmExtraPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"

//Optical processes
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

//Decays
#include "G4Decay.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecayBase.hh"
#include "G4NuclideTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSPhysicsList::SAXSPhysicsList():
  G4VUserPhysicsList(),
  fUseMIFlag(true),
  fPMessenger(0)
{
  G4cout << "### PhysicsList instantiated ###" << G4endl; 
  
  G4LossTableManager::Instance();
  
  //set default cuts value
  defaultCutValue = 0.1*mm;
  
  //define the messenger
  fPMessenger = new SAXSPhysicsListMessenger(this);
  
  //set verbosity
  SetVerboseLevel(1);
  
  //particle list
  fParticleList = new G4DecayPhysics(verboseLevel);
  
  //EM Physics
  fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSPhysicsList::~SAXSPhysicsList() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSPhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSPhysicsList::ConstructProcess()
{        
  //transportation
  AddTransportation();
  
  //EM physics
  fEmPhysicsList->ConstructProcess();

  //Atomic deexcitation
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(true);                 //Activate deexcitation processes and fluorescence
  de->SetAuger(true);            //Activate Auger effect if deexcitation is activated
  de->SetAugerCascade(true); //Activate Auger Cascade if deexcitation is activated
  de->SetPIXE(true);                   //Activate Particle Induced X-Ray Emission (PIXE)
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSPhysicsList::SelectPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "### PhysicsList::SelectPhysicsList: <" << name << "> ###" << G4endl;
  }

  if (name == "standard") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
    G4cout << "### selected Standard PhysicsList ###" << G4endl;
  } else if (name == "standard_option4") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);
    G4cout << "### selected Standard_option4 PhysicsList ###" << G4endl;
  } else if (name == "livermore") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);
    G4cout << "### selected Livermore PhysicsList ###" << G4endl;
  } else if (name == "penelope") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);
    G4cout << "### selected Penelope PhysicsList ###" << G4endl;
  } else if (name == "penelopeMI") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysicsMI(verboseLevel,"G4EmPenelopeMI",fUseMIFlag);
    G4cout << "### selected Penelope PhysicsList with MI effects ###" << G4endl;
  } else if (name == "LowEP") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics(1,name);
    G4cout << "### selected LowEP PhysicsList ###" << G4endl;            
  } else {
    G4cout << "### PhysicsList::SelectPhysicsList: <" << name 
           << ">"<< " is not defined ###" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSPhysicsList::SetCuts()
{
  //set the default cuts value for all particle types  
  SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSPhysicsList::SetDefaultCutsValue(G4double value)
{
  //define a new default cuts value
  defaultCutValue = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

