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
/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGenerator::~PrimaryGenerator()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
  //vertex A uniform on a cylinder
  //
  const G4double r = 2*mm;
  const G4double zmax = 8*mm;
  //
  G4double alpha = twopi*G4UniformRand();     //alpha uniform in (0, 2*pi)
  G4double ux = std::cos(alpha);
  G4double uy = std::sin(alpha);
  G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)
  G4ThreeVector positionA(r*ux,r*uy,z);
  G4double timeA = 0*s;
  // 
  G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);
  
  //particle 1 at vertex A
  //
  G4ParticleDefinition* particleDefinition
           = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
  particle1->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
  particle1->SetKineticEnergy(1*MeV);
  //
  vertexA->SetPrimary(particle1);
  event->AddPrimaryVertex(vertexA);

  //vertex (B) symetric to vertex A
  //
  alpha += pi;
  ux = std::cos(alpha);
  uy = std::sin(alpha);
  G4ThreeVector positionB(r*ux,r*uy,z);
  G4double timeB = 1*s;
  // 
  G4PrimaryVertex* vertexB = new G4PrimaryVertex(positionB, timeB);
  
  //particles 2 and 3 at vertex B
  //
  const G4double dalpha = 10*deg;
  ux = std::cos(alpha + dalpha);
  uy = std::sin(alpha + dalpha);        
  G4PrimaryParticle* particle2 = new G4PrimaryParticle(particleDefinition);
  particle2->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
  particle2->SetKineticEnergy(1*keV);
  //
  ux = std::cos(alpha - dalpha);
  uy = std::sin(alpha - dalpha);        
  G4PrimaryParticle* particle3 = new G4PrimaryParticle(particleDefinition);
  particle3->SetMomentumDirection(G4ThreeVector(ux,uy,0));    
  particle3->SetKineticEnergy(1*GeV);
  //
  vertexB->SetPrimary(particle2);
  vertexB->SetPrimary(particle3);  
  event->AddPrimaryVertex(vertexB);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
