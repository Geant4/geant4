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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option2.hh"
<<<<<<< HEAD
=======
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4DNAGenericIonsManager.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
: G4VModularPhysicsList()
{
  SetDefaultCutValue(1.0*nm);
  SetVerboseLevel(1);
  
  //RegisterPhysics(new G4EmDNAPhysics());
  RegisterPhysics(new G4EmDNAPhysics_option2());
  
<<<<<<< HEAD
  G4ProductionCutsTable::GetProductionCutsTable()->
      SetEnergyRange(1*eV, 1*GeV);
=======
  // electromagnetic physics list
  //
  fEmPhysicsList->ConstructProcess();
      
  // tracking cut
  //
  AddTrackingCut();

  // maximum step size
  //
  AddMaxStepSize();

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
PhysicsList::~PhysicsList()
=======
void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>-1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "dna") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics();
         
  } else if (name == "dna_opt1") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option1();
         
  } else if (name == "dna_opt2") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option2();
         
  } else if (name == "dna_opt3") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option3();
         
  } else if (name == "dna_opt4") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option4();
         
  } else if (name == "dna_opt5") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option5();
         
  } else if (name == "dna_opt6") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option6();
         
  } else if (name == "liv") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
         
  } else if (name == "pene") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
         
  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddTrackingCut()
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddMaxStepSize()
{

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") 
    {
      ph->RegisterProcess(new G4StepLimiter(), particle); 
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
