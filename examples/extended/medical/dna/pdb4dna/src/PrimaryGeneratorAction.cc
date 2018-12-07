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
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fpParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic

  G4ParticleDefinition* particle
  = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fpParticleGun->SetParticleDefinition(particle);
  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fpParticleGun->SetParticleEnergy(0.1*MeV);
  fpParticleGun->SetParticlePosition(G4ThreeVector(0.*nm,0.*nm,0.*nm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  G4double boundXHalfLength = 0;
  G4double boundYHalfLength = 0;
  G4double boundZHalfLength = 0;

  G4VPhysicalVolume* boundPV
  = G4PhysicalVolumeStore::GetInstance()->GetVolume("boundingPV");

  G4ThreeVector boundPos;
  if ( boundPV )
  {
    boundPos = boundPV->GetTranslation();
  }

  G4LogicalVolume* boundLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingLV");

  G4Box* boundBox = 0;
  if ( boundLV )
  {
    boundBox = dynamic_cast< G4Box*>(boundLV->GetSolid());
  }

  if ( boundBox )
  {
    boundXHalfLength = boundBox->GetXHalfLength(); 
    boundYHalfLength = boundBox->GetYHalfLength(); 
    boundZHalfLength = boundBox->GetZHalfLength();

    // Set gun position
    // Select a starting position on a sphere including the target volume
    //
    G4double radius = std::sqrt(boundXHalfLength*boundXHalfLength+
        boundYHalfLength*boundYHalfLength+
        boundZHalfLength*boundZHalfLength);
    G4double cosTheta = 2*G4UniformRand()-1;
    G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
    G4double phi = twopi*G4UniformRand();

    G4ThreeVector positionStart(boundPos.x()+radius*sinTheta*std::cos(phi),
        boundPos.y()+radius*sinTheta*std::sin(phi),
        boundPos.z()+radius*cosTheta);

    fpParticleGun->SetParticlePosition(positionStart);

    // Set gun direction
    // To compute the direction, select a point inside the target volume
    //
    G4ThreeVector positionDir(
        boundPos.x()+boundXHalfLength*(2*G4UniformRand()-1),
        boundPos.y()+boundYHalfLength*(2*G4UniformRand()-1),
        boundPos.z()+boundZHalfLength*(2*G4UniformRand()-1));

    fpParticleGun->SetParticleMomentumDirection(
        (positionDir-positionStart).unit());
  }
  else
  {
    G4cerr << "Bounding volume not found." << G4endl;
    G4cerr << "Default particle kinematic used" << G4endl;
  } 

  fpParticleGun->GeneratePrimaryVertex(anEvent);
}
