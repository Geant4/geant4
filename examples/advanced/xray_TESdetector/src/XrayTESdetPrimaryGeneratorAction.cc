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
/// \file XrayTESdetPrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetPrimaryGeneratorAction.hh"
#include "XrayTESdetDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayTESdetPrimaryGeneratorAction::XrayTESdetPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  // setting particle type
   fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="proton"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayTESdetPrimaryGeneratorAction::~XrayTESdetPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double ax;
  G4double ay;
  G4double az;

  G4double bx;
  G4double by;
  G4double bz;

  G4double vx;
  G4double vy;
  G4double vz;

  G4double energy = 0;
  G4double u;
  G4double u1;

  G4double theta0;
  G4double phi0;
  G4double theta1;
  G4double phi1;

  G4double ene;
  G4double flux;
  G4double r;

  // Radius to shot particles
	r = 27*cm;

  // Proton flux @L2 (Pamela2009), in p/cm2/s/sr [10 MeV - 10 GeV] = 0.407
  G4double phi = 0.379;
  G4double tr = 0.938;
  G4double T;

  while (true)
  {
    u = G4UniformRand();
    u1 = 0.268*G4UniformRand();      //max of the flux, in particles/cm2-s-GeV-sr
    ene = 10.*MeV + u*100000.*MeV;
    T = ene*0.001;
    flux = 1.9*std::pow(((T + phi)*(T + phi + 2*tr)), -1.39)*T*(T + 2*tr)/((1 + 0.4866*std::pow(((T + phi)*(T + phi + 2*tr)),-1.255))*(T + phi)*(T + phi + 2*tr));
    if (u1 <= flux)
    {
      energy = ene;
      break;
    }
  }

  theta0 = std::acos(1. - 2.*G4UniformRand());
  theta1 = std::acos(1. - 2.*G4UniformRand());
  phi0 = twopi*G4UniformRand();
  phi1 = twopi*G4UniformRand();
  ax = r*std::sin(theta0)*std::cos(phi0);
  ay = r*std::sin(theta0)*std::sin(phi0);
  az = r*std::cos(theta0) + 80*mm;

  bx = r*std::sin(theta1)*std::cos(phi1);
  by = r*std::sin(theta1)*std::sin(phi1);
  bz = r*std::cos(theta1) + 80*mm;
  vx = bx - ax;
  vy = by - ay;
  vz = bz - az;

  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(G4ThreeVector(ax, ay, az));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
