// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em6PrimaryGeneratorAction.cc,v 1.3 2000-01-28 04:26:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em6PrimaryGeneratorAction.hh"

#include "Em6DetectorConstruction.hh"
#include "Em6PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 G4String Em6PrimaryGeneratorAction::thePrimaryParticleName="e-" ; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6PrimaryGeneratorAction::Em6PrimaryGeneratorAction(
                                            Em6DetectorConstruction* Em6DC)
:Em6Detector(Em6DC)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  
  //create a messenger for this class
  gunMessenger = new Em6PrimaryGeneratorMessenger(this);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6PrimaryGeneratorAction::~Em6PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;  
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6PrimaryGeneratorAction::SetDefaultKinematic()
{    
  // default particle kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  
  thePrimaryParticleName = particle->GetParticleName();

  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(30.*MeV);
  //  G4double x0 = -0.5*(Em6Detector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(-0.1*cm, 0*cm, 0*cm));  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  thePrimaryParticleName = particleGun->GetParticleDefinition()->
                                                GetParticleName();
   
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String Em6PrimaryGeneratorAction::GetPrimaryName()
{
   return thePrimaryParticleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

