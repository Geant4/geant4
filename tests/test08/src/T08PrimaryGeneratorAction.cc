// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PrimaryGeneratorAction.cc,v 1.1 1999-01-08 16:35:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// $ Id: $
//
// From exampleEmPhys2/MyPrimaryGeneratorAction.cc,v 1.1 1998/02/06  maire 

#include "T08PrimaryGeneratorAction.hh"

#include "T08DetectorConstruction.hh"
#include "T08PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

T08PrimaryGeneratorAction::T08PrimaryGeneratorAction(T08DetectorConstruction* myDC)
:myDetector(myDC),rndmFlag("off")
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new T08PrimaryGeneratorMessenger(this);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="proton");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(3.0*GeV);
  G4double position = -0.5*(myDetector->GetWorldFullLength());
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));

}

T08PrimaryGeneratorAction::~T08PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void T08PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 
  // The Target is a box placed at (0,0,0)
  //
  G4double x0 = -0.5*(myDetector->GetTargetFullLength());
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (myDetector->GetTargetFullLength())*(G4UniformRand()-0.5);
      z0 = (myDetector->GetTargetFullLength())*(G4UniformRand()-0.5);
     } 
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  particleGun->GeneratePrimaryVertex(anEvent);
}


