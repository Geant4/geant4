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
/// \file field/field06/src/F06PrimaryGeneratorAction.cc
/// \brief Implementation of the F06PrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F06PrimaryGeneratorAction.hh"

#include "F06DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06PrimaryGeneratorAction::F06PrimaryGeneratorAction(void)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  G4double mass_proton = particle->GetPDGMass();

  G4double muN = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/
                                               (mass_proton/CLHEP::c_squared);

  G4cout << " *** Neutron *** " << G4endl;

  particle = particleTable->FindParticle("neutron");
  G4double mass_neutron = particle->GetPDGMass();

  G4double magneticMoment = particle->GetPDGMagneticMoment();
  G4cout << " magneticMoment: " << magneticMoment/muN << G4endl;

  // g_factor for spin 1/2

  G4double g_factor = 2 * magneticMoment/muN;
  G4cout << " g_factor: " << g_factor << G4endl;

  G4double charge = particle->GetPDGCharge();
  G4cout << " charge: " << charge << G4endl;

  G4double anomaly = (g_factor - 2.)/2.;
  G4cout << " anomaly: " << anomaly << G4endl;

  anomaly = (g_factor * (mass_neutron/mass_proton) - 2.)/2.;
  G4cout << " corrected anomaly: " << anomaly << G4endl;

  muN = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/(mass_neutron/CLHEP::c_squared);
  g_factor = 2 * magneticMoment/muN;

  anomaly = (g_factor - 2.)/2.;
  G4cout << " *** anomaly: " << anomaly << G4endl;

  G4cout << " *** MuonPlus *** " << G4endl;

  particle = particleTable->FindParticle("mu+");
  G4double mass_muon = particle->GetPDGMass();

  G4double muB = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/
                                                 (mass_muon/CLHEP::c_squared);

  magneticMoment = particle->GetPDGMagneticMoment();
  G4cout << " magneticMoment: " << magneticMoment/muB << G4endl;

  // g_factor for spin 1/2

  g_factor = magneticMoment/muB;
  G4cout << " g_factor: " << g_factor << G4endl;

  charge = particle->GetPDGCharge();
  G4cout << " charge: " << charge << G4endl;

  anomaly = (g_factor - 2.)/2.;
  G4cout << " anomaly: " << anomaly << G4endl;

//  anomaly = particle->CalculateAnomaly();
//  G4cout << " *** anomaly: " << anomaly << G4endl;

  particle = particleTable->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06PrimaryGeneratorAction::~F06PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F06PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //

  fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
  fParticleGun->SetParticlePolarization(G4ThreeVector(0,1,0));

  G4double particleEnergy = G4UniformRand()*1e-7*eV;
  fParticleGun->SetParticleEnergy(particleEnergy);

  G4double theta = 2*pi*G4UniformRand();
  G4double phi = std::acos(1-2*G4UniformRand());
  if (phi > pi/2 && phi < pi) phi = pi-phi;

  G4double z = std::sin(phi)*std::cos(theta);
  G4double x = std::sin(phi)*std::sin(theta);
  G4double y = std::cos(phi);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
