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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SingleParticleSource.hh"
#include "G4ParticleTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction *pPrimaryGenerator)
    : G4UImessenger(), fpPrimaryGenerator(pPrimaryGenerator) {
  fparticle = std::make_unique<G4UIcmdWithAString>("/UHDR/source/particle", this);
  fparticle->SetGuidance("Add time structure.");
  fparticle->SetGuidance("e-, proton, alpha");

  fenergy = std::make_unique<G4UIcmdWithADoubleAndUnit>("/UHDR/source/energy", this);
  fenergy->SetGuidance("Sets a monocromatic energy (same as  gps/energy)");
  fenergy->SetParameterName("monoenergy", false, false);
  fenergy->SetDefaultUnit("keV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
  auto particleGun = fpPrimaryGenerator->GetSPGun();
  if (command == fparticle.get()) {
    G4ExceptionDescription ed;
    G4ParticleDefinition *pd = G4ParticleTable::GetParticleTable()->FindParticle(newValue);
    if (pd != nullptr) {
      particleGun->SetParticleDefinition(pd);
    } else {
      ed << "Particle [" << newValue << "] is not found.";
      command->CommandFailed(ed);
    }
  } else if (command == fenergy.get()) {
    auto pEnergyDis = particleGun->GetEneDist();
    pEnergyDis->SetMonoEnergy(fenergy->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......