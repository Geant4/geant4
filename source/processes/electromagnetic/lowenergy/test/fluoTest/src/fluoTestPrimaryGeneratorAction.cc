//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



#include "fluoTestPrimaryGeneratorAction.hh"

#include "fluoTestDetectorConstruction.hh"
#include "fluoTestPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestPrimaryGeneratorAction::fluoTestPrimaryGeneratorAction(
                                               fluoTestDetectorConstruction* DC)
  :Detector(DC),rndmFlag("off"),randomizePrimary("off"),sigmaAngle(0),sigmaMomentum(0)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new fluoTestPrimaryGeneratorMessenger(this);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                   = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(50.*MeV);
  
  G4double position = -0.5*(Detector->GetWorldSizeX());
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestPrimaryGeneratorAction::~fluoTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the beginning of each event
  
  G4ParticleDefinition* particle;

  G4double x0 = -0.5*(Detector->GetWorldSizeX());
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (Detector->GetSampleSizeYZ())*(G4UniformRand()-0.5);
      z0 = (Detector->GetSampleSizeYZ())*(G4UniformRand()-0.5);
     particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
     }
  //particleGun->GeneratePrimaryVertex(anEvent);
 
 if(randomizePrimary =="on")
  {
    G4double pp = momentum + (G4UniformRand()-0.5)*sigmaMomentum;
    G4double mass = particle->GetPDGMass();
    G4double  Ekin = sqrt(pp*pp+mass*mass)-mass;
    particleGun->SetParticleEnergy(Ekin);
    G4double angle = (G4UniformRand()-0.5)*sigmaAngle;
    particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(angle),0.,cos(angle)));
 
  }
particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestPrimaryGeneratorAction::SetSigmaMomentum(G4double val)
{
  sigmaMomentum = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void fluoTestPrimaryGeneratorAction::SetSigmaAngle(G4double val)
{
  sigmaAngle = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....










