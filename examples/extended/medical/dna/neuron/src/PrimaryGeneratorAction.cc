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
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $ID$
/// \file PrimaryGeneratorAction.cc 
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
//
#include "CommandLineParser.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

using namespace G4DNAPARSER ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction():
G4VUserPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fpParticleGun  = new G4ParticleGun(n_particle);
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("proton");
  fpParticleGun->SetParticleDefinition(particle);
  // default gun parameters
  fpParticleGun->SetParticleEnergy(10.*MeV);   
  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fpParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.*um));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   
  // Initial kinetic energy of particles beam as Gaussion distribution! 
  
  // G4double Ekin = 10.*MeV;
  // G4double deltaEkin = 9.0955e-5*MeV;   
  // fpParticleGun->SetParticleEnergy(G4RandGauss::shoot(Ekin,deltaEkin)); 
   
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore

  // We included three options for particles direction:
  
  G4double mediumRadius     = 0.;
  G4double boundingXHalfLength = 0.;
  G4double boundingYHalfLength = 0.;
  G4double boundingZHalfLength = 0.; 
  G4LogicalVolume* mediumLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume("Medium");  
  G4LogicalVolume* boundingLV
  = G4LogicalVolumeStore::GetInstance()->GetVolume("BoundingSlice");  
  G4Orb* mediumSphere  = 0;
  G4Box* boundingSlice = 0;  
  if ( mediumLV && boundingLV)
  {
   mediumSphere = dynamic_cast< G4Orb*>(mediumLV->GetSolid());
   boundingSlice = dynamic_cast< G4Box*>(boundingLV->GetSolid());
  } 
  if ( mediumSphere && boundingSlice)
  {
    mediumRadius   = mediumSphere->GetRadius();
    boundingXHalfLength = boundingSlice->GetXHalfLength(); 
    boundingYHalfLength = boundingSlice->GetYHalfLength(); 
    boundingZHalfLength = boundingSlice->GetZHalfLength();
 
  /// a) Partilces directed to "square" on the XY plane of bounding slice 
  ///   (or YZ, XZ) 
 
    if ( CommandLineParser::GetParser()->GetCommandIfActive("-sXY"))
     {
   //G4cerr << " Initial beam position uniformly spread on a square! " 
   //      << G4endl;
    // INITIAL BEAM DIVERGENCE
    fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); 
    // // INITIAL BEAM POSITION
    fpParticleGun->SetParticlePosition(G4ThreeVector(
    CLHEP::RandFlat::shoot(-boundingXHalfLength,boundingXHalfLength),
    CLHEP::RandFlat::shoot(-boundingYHalfLength,boundingYHalfLength),
    -mediumRadius));
   // Surface area of a square 
    // fGunArea = 4*boundingXHalfLength*boundingYHalfLength ;
    //G4cerr << " Particle Fluence Area on a Square (um2) =  " 
    //       <<fGunArea / (um*um)<< G4endl;
     }

  /// b) Partilces directed to "disk" on the XY plane of 
  ///    bounding slice (or YZ, XZ) 

    else if ( CommandLineParser::GetParser()->GetCommandIfActive("-dXY"))
     {
    //G4cerr << " Initial beam position uniformly spread on a disk! " 
    //       << G4endl;
    fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.)); 
    G4double x0,y0,z0;
    x0 = 100.*mm;
    y0 = 100.*mm;
    z0 = -mediumRadius; // mediumRadius;
    while (! (std::sqrt(x0*x0+y0*y0)<= mediumRadius) )
    {
     x0 = CLHEP::RandFlat::shoot(-mediumRadius,mediumRadius);
     y0 = CLHEP::RandFlat::shoot(-mediumRadius,mediumRadius); 
    }  
    fpParticleGun->SetParticlePosition(G4ThreeVector(x0, y0,z0));     
   // Surface area of a disk 
   // fGunArea = pi*mediumRadius*mediumRadius ;
   // G4cerr << " Particle Fluence Area on a Disk (um2) =  " 
   //        <<fGunArea / (um*um)<< G4endl;

     } 
 
  /// c) Partilces directed towards the bounding slice (default option!) 
   // Select a starting position on a sphere including the 
   // target volume and neuron morphology
   else 
   {
   //G4cerr << " Initial beam position uniformly spread on a sphere! "<< G4endl;
    G4double cosTheta = 2.*G4UniformRand()-1;
    G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
    G4double phi = twopi*G4UniformRand();
    G4ThreeVector positionStart(mediumRadius*sinTheta*std::cos(phi),
        mediumRadius*sinTheta*std::sin(phi),
        mediumRadius*cosTheta);  
    fpParticleGun->SetParticlePosition(positionStart); 
    // To compute the direction, select a point inside the target volume
    G4ThreeVector positionDir(
        boundingXHalfLength*(2.*G4UniformRand()-1),
        boundingYHalfLength*(2.*G4UniformRand()-1),
        boundingZHalfLength*(2.*G4UniformRand()-1));
    fpParticleGun->SetParticleMomentumDirection(
        (positionDir-positionStart).unit());
   // Surface area of sphere 
    fGunArea = 4.*pi*mediumRadius*mediumRadius ;
   //  G4cerr << " Particle Fluence Area on sphere (um2) =  " 
   //         <<fGunArea / (um*um)<< G4endl; 
 }
  
  }
  else
  {
    G4cerr << "Bounding slice volume not found!" << G4endl;
    G4cerr << "Default particle kinematic used" << G4endl;
  }  

  fpParticleGun->GeneratePrimaryVertex(anEvent);
}
