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
/// \file electromagnetic/TestEm17/src/PhysListEmStandard.cc
/// \brief Implementation of the PhysListEmStandard class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 

#include "PhysListEmStandard.hh"

#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Proton.hh"

#include "G4ProcessManager.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetMinEnergy(100*eV);  
  param->SetMaxEnergy(1000*PeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  // mu+
  G4ParticleDefinition* particle = G4MuonPlus::MuonPlus();
  G4ProcessManager* pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4MuIonisation(),     -1, 2, 2);
  pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
  pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);

  // mu-
  particle = G4MuonMinus::MuonMinus();
  pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4MuIonisation(),     -1, 2, 2);
  pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
  pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);

  // pi+
  particle = G4PionPlus::PionPlus();
  pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4hIonisation(),     -1, 2, 2);
  pmanager->AddProcess(new G4hBremsstrahlung(), -1, 3, 3);
  pmanager->AddProcess(new G4hPairProduction(), -1, 4, 4);

  // pi-
  particle = G4PionMinus::PionMinus();
  pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4hIonisation(),     -1, 2, 2);
  pmanager->AddProcess(new G4hBremsstrahlung(), -1, 3, 3);
  pmanager->AddProcess(new G4hPairProduction(), -1, 4, 4);

  // proton
  particle = G4Proton::Proton();
  pmanager = particle->GetProcessManager();    
  pmanager->AddProcess(new G4hIonisation(),     -1, 2, 2);
  pmanager->AddProcess(new G4hBremsstrahlung(), -1, 3, 3);
  pmanager->AddProcess(new G4hPairProduction(), -1, 4, 4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

