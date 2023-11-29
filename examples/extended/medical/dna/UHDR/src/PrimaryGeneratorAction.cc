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

#include "PrimaryGeneratorAction.hh"
#include "ChemistryWorld.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SingleParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Electron.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction *pDet)
    : G4VUserPrimaryGeneratorAction(), fpDetector(pDet) {
  fpMessenger = std::make_unique<PrimaryGeneratorMessenger>(this);
  fParticleGun = std::make_unique<G4SingleParticleSource>();
  G4ParticleDefinition *particle = G4Electron::Definition();
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetNumberOfParticles(1000000); // by user

  auto pPosDist = fParticleGun->GetPosDist();
  pPosDist->SetPosDisType("Plane");
  pPosDist->SetPosDisShape("Square");
  auto faceSiez = fpDetector->GetChemistryWorld()
      ->GetChemistryBoundary()
      ->halfSideLengthInY();
  pPosDist->SetCentreCoords(G4ThreeVector(0, 0, -faceSiez));
  pPosDist->SetHalfX(faceSiez);
  pPosDist->SetHalfY(faceSiez);

  auto pAngleDist = fParticleGun->GetAngDist();
  pAngleDist->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  auto pEnergyDis = fParticleGun->GetEneDist();
  pEnergyDis->SetMonoEnergy(0.9999 * MeV);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  G4ParticleDefinition *particle = fParticleGun->GetParticleDefinition();
  auto NumberOfParticlesToBeGenerated = fParticleGun->GetNumberOfParticles();
  auto pPosDist = fParticleGun->GetPosDist();
  auto pAngleDist = fParticleGun->GetAngDist();
  auto pEnDist = fParticleGun->GetEneDist();
  auto rnd = fParticleGun->GetBiasRndm();
  auto charge = particle->GetPDGCharge();

  for (G4int i = 0; i < NumberOfParticlesToBeGenerated; i++) {
    auto angle = pAngleDist->GenerateOne();
    auto energy = pEnDist->GenerateOne(particle);
    auto pos = pPosDist->GenerateOne();
    auto mass = particle->GetPDGMass();
    auto p = new G4PrimaryParticle(particle);
    auto vertex = new G4PrimaryVertex(pos, 0);
    p->SetKineticEnergy(energy);
    p->SetMass(mass);
    p->SetMomentumDirection(angle);
    p->SetCharge(charge);
    p->SetPolarization(0., 0., 0.);

    G4double weight = pEnDist->GetWeight() * rnd->GetBiasWeight();
    p->SetWeight(weight);
    vertex->SetPrimary(p);
    anEvent->AddPrimaryVertex(vertex);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
