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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "PhysicsList.hh"

#include "G4ProcessManager.hh"

#include "G4DecayPhysics.hh"
#include "G4EmDNABuilder.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4GenericIon.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclideTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
{
  SetDefaultCutValue(1.0 * micrometer);
  SetVerboseLevel(1);

  fEmPhysics = "DNA_Opt2";
  fEmPhysicsList = new G4EmDNAPhysics_option2();
  fDecayPhysicsList = new G4DecayPhysics();

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV, 1 * GeV);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetMinEnergy(100 * eV);
  param->SetMaxEnergy(1 * GeV);

  // Limits in G4NuclideTable
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1 * picosecond);
  G4NuclideTable::GetInstance()->SetLevelTolerance(1.0 * eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fEmPhysicsList;
  delete fDecayPhysicsList;
  delete fRadDecayPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fEmPhysicsList->ConstructParticle();
  fDecayPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  fEmPhysicsList->ConstructProcess();

  fDecayPhysicsList->ConstructProcess();

  if (nullptr != fRadDecayPhysicsList) fRadDecayPhysicsList->ConstructProcess();

  if (fIsTrackingCutSet) TrackingCut();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysics(const G4String& name)
{
  if (name == fEmPhysics) {
    return;
  }

  G4cout << "### PhysicsList::AddPhysics Warning: Physics List <" << name << "> is requested"
         << G4endl;

  fEmPhysics = name;

  if (name == "emstandard_opt0") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
  }
  else if (name == "emstandard_opt3") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
  }
  else if (name == "emstandard_opt4") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
  }
  else if (name == "raddecay") {
    if (nullptr == fRadDecayPhysicsList) fRadDecayPhysicsList = new G4RadioactiveDecayPhysics();
  }
  else if (name == "emlivermore") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
  }
  else if (name == "empenelope") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
  }
  else if (name == "DNA_Opt0") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics();
  }
  else if (name == "DNA_Opt1") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option1();
  }
  else if (name == "DNA_Opt2") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option2();
  }
  else if (name == "DNA_Opt3") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option3();
  }
  else if (name == "DNA_Opt4") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option4();
  }
  else if (name == "DNA_Opt5") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option5();
  }
  else if (name == "DNA_Opt6") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option6();
  }
  else if (name == "DNA_Opt7") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option7();
  }
  else if (name == "DNA_Opt8") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmDNAPhysics_option8();
  }
  else {
    G4cout << "### PhysicsList::AddPhysics Warning: Physics List <" << name
           << "> is does not exist - the command ignored" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::TrackingCut()
{
  auto particle = G4GenericIon::GenericIon();  // DNA heavy ions
  auto particleName = particle->GetParticleName();
  auto capture = G4EmDNABuilder::FindOrBuildCapture(0.5 * CLHEP::MeV, particle);
  capture->AddRegion("World");
  capture->SetKinEnergyLimit(0.5 * CLHEP::MeV);  // 0.5 MeV/u
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetTrackingCut(G4bool isCut)
{
  fIsTrackingCutSet = isCut;
}
