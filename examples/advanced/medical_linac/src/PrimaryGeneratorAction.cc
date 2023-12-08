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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include <G4ParticleTable.hh>
#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleGun.hh>
#include <G4GeneralParticleSource.hh>
#include <Randomize.hh>
#include <G4RunManager.hh>
#include <G4PhysicalConstants.hh>

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction()
: fGun(nullptr),
  gunRadius(0), gunMeanEnergy(0), gunStdEnergy(0),
  energy(0),
  particle(nullptr)
{
  priGenMessenger = new PrimaryGeneratorMessenger(this);  
  G4int n_particle = 1;
  fGun = new G4ParticleGun(n_particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete priGenMessenger;
    delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  GenerateFromRandom();
  
  fGun -> SetParticlePosition(position);
  fGun -> SetParticleMomentumDirection((G4ParticleMomentum)direction);
  fGun -> SetParticleEnergy(energy * MeV);
  fGun -> SetParticleDefinition(particle);
  // create vertex with previous specifications
  fGun -> GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::GenerateFromRandom()
{
  DetectorConstruction* detector = (DetectorConstruction*)
    G4RunManager::GetRunManager() -> GetUserDetectorConstruction();
  
  G4double zPos = -(detector -> GetAccOriginPosition()) - 5.*mm; // 5mm before the target disk
  G4double alpha, rho, phi, sinTheta, cosTheta;
  
  sinTheta = G4RandGauss::shoot(0., 0.003);
  cosTheta = sqrt(1-sinTheta*sinTheta);
  phi = twopi * G4UniformRand();
  rho = gunRadius * G4UniformRand();
  alpha = twopi * G4UniformRand();
  
  particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  direction = {sinTheta * cos(phi), sinTheta * sin(phi), cosTheta};
  position = {rho * sin(alpha), rho * cos(alpha), zPos};
  energy = G4RandGauss::shoot(gunMeanEnergy, gunStdEnergy);
}

void PrimaryGeneratorAction::SetGunRadius(G4double value)
{
  gunRadius = value;
}

void PrimaryGeneratorAction::SetGunMeanEnergy(G4double value)
{
  gunMeanEnergy = value;
}

void PrimaryGeneratorAction::SetGunStdEnergy(G4double value)
{
  gunStdEnergy = value;
}
