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
// $Id$
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Sc01PrimaryGeneratorAction.hh"
#include "Sc01DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <iostream>
#include <fstream>

Sc01PrimaryGeneratorAction::
Sc01PrimaryGeneratorAction(Sc01DetectorConstruction*det)
  :
  fDetector(det)
{
  particleGun = new G4ParticleGun();
}

Sc01PrimaryGeneratorAction::~Sc01PrimaryGeneratorAction()
{
  delete particleGun;
}

void Sc01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* aParticleDefinition 
    = particleTable->FindParticle(particleName="opticalphoton");

  // create a new vertex in a random position 

  //  G4ThreeVector VertexPosition (GetRandomPosition());
  G4ThreeVector VertexPosition (G4ThreeVector(0,0,0));
  G4PrimaryVertex* aVertex = 
    new G4PrimaryVertex( VertexPosition, 0);

  //  G4int NumberOfParticlesToBeGenerated = 10000;
  G4int NumberOfParticlesToBeGenerated = 1;
  // G4cout << "  A " << NumberOfParticlesToBeGenerated 
  //          << " optical photons vertex has been generated at " 
  //          << VertexPosition << G4endl;
  // create new primaries and set them to the vertex
  for( int i=0; i<NumberOfParticlesToBeGenerated; i++ )
  {
    //	G4ThreeVector m = GetRandomDirection();
    G4ThreeVector rm(0.,1.,0.);
    G4PrimaryParticle* aPrimaryParticle =
      new G4PrimaryParticle(aParticleDefinition, rm.x(), rm.y(), rm.z());
    aPrimaryParticle->SetMass (0);
    G4ThreeVector p = GetRandomPolarization ( rm );
    aPrimaryParticle->SetPolarization(p.x(),p.y(),p.z());
    aVertex->SetPrimary( aPrimaryParticle );
  }
  anEvent->AddPrimaryVertex( aVertex );
}

G4ThreeVector Sc01PrimaryGeneratorAction::GetRandomDirection() {

  G4ThreeVector retval;

  G4double CosTheta;
  G4double SinTheta;

  G4double Phi;
  G4double SinPhi;
  G4double CosPhi;

  G4double rand;

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = std::sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = std::sin (Phi);
  CosPhi = std::cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}

G4ThreeVector Sc01PrimaryGeneratorAction::GetRandomPosition() 
{
  G4double a = fDetector->GetHallSize();

  G4double x = ( G4UniformRand()*2 - 1 )*a;
  G4double y = ( G4UniformRand()*2 - 1 )*a;
  G4double z = ( G4UniformRand()*2 - 1 )*a;

  G4ThreeVector retval (x, y, z);

  return retval;
}

G4ThreeVector Sc01PrimaryGeneratorAction::GetRandomPolarization(G4ThreeVector Direction) {
  G4ThreeVector Polarization = Direction.orthogonal();
  G4ThreeVector retval = Polarization.unit();
  return retval;
}






