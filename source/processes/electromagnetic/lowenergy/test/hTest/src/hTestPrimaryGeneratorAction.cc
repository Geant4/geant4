// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: hTestPrimaryGeneratorAction.cc,v 1.4 2001-03-12 09:17:15 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestPrimaryGeneratorAction.hh"

#include "hTestDetectorConstruction.hh"
#include "hTestPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 G4String hTestPrimaryGeneratorAction::thePrimaryParticleName="e-" ; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::hTestPrimaryGeneratorAction(
                                            hTestDetectorConstruction* hTestDC)
:hTestDetector(hTestDC)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  
  //create a messenger for this class
  gunMessenger = new hTestPrimaryGeneratorMessenger(this);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorAction::~hTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;  
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::SetDefaultKinematic()
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
  //  G4double x0 = -0.5*(hTestDetector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm));  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 

  thePrimaryParticleName = particleGun->GetParticleDefinition()->
                                                GetParticleName();
   
  G4double energy;
  //  G4double sigmaE = 0.9*MeV ;
  G4double sigmaE = -0.9*MeV ;

  if(0.0 < sigmaE) {
    energy += G4RandGauss::shoot(175.0*MeV,sigmaE);
    if(0.0 > energy) energy = 0.0;
    particleGun->SetParticleEnergy(energy);
  }  

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String hTestPrimaryGeneratorAction::GetPrimaryName()
{
   return thePrimaryParticleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

