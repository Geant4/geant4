// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2PrimaryGeneratorAction.cc,v 1.2 1999-12-15 14:49:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2PrimaryGeneratorAction.hh"

#include "Em2DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2PrimaryGeneratorAction::Em2PrimaryGeneratorAction(
                                               Em2DetectorConstruction* det)
:Em2Detector(det)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  G4ParticleDefinition* particle
                 = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(5.*GeV);
  G4double position = -0.5*(Em2Detector->GetfullLength());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2PrimaryGeneratorAction::~Em2PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  G4double position = -0.5*(Em2Detector->GetfullLength());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));     
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

