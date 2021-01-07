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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction1.cc
/// \brief Implementation of the PrimaryGeneratorAction1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction1.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction1::PrimaryGeneratorAction1(G4ParticleGun* gun)
: fParticleGun(gun)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction1::~PrimaryGeneratorAction1()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction1::GeneratePrimaries(G4Event* anEvent)
{
  const G4double r = 2*mm;
  const G4double zmax = 8*mm;   
  
  //vertex 1 uniform on cylinder
  //
  G4double alpha = twopi*G4UniformRand();  //alpha uniform in (0, 2*pi)
  G4double ux = std::cos(alpha);
  G4double uy = std::sin(alpha);
  G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)
        
  fParticleGun->SetParticlePosition(G4ThreeVector(r*ux,r*uy,z));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  fParticleGun->SetParticleEnergy(1*MeV);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  //vertex 2 at opposite
  //
  alpha += pi;
  ux = std::cos(alpha);
  uy = std::sin(alpha);        
  fParticleGun->SetParticlePosition(G4ThreeVector(r*ux,r*uy,z));
  
  //particle 2 at vertex 2
  //
  const G4double dalpha = 10*deg;
  ux = std::cos(alpha + dalpha);
  uy = std::sin(alpha + dalpha);        
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  fParticleGun->SetParticleEnergy(1*keV);
  fParticleGun->GeneratePrimaryVertex(anEvent);
   
  //particle 3 at vertex 3 (same as vertex 2)
  //
  ux = std::cos(alpha - dalpha);
  uy = std::sin(alpha - dalpha);        
  
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  fParticleGun->SetParticleEnergy(1*GeV);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  // randomize time zero of anEvent
  //
  G4double tmin = 0*s, tmax = 10*s;
  G4double t0 = tmin + (tmax - tmin)*G4UniformRand();
  G4PrimaryVertex* aVertex = anEvent->GetPrimaryVertex();
  while (aVertex) {
    aVertex->SetT0(t0);
    aVertex = aVertex->GetNext();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
