//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



#include "XrayFluoPrimaryGeneratorAction.hh"

#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPrimaryGeneratorAction::XrayFluoPrimaryGeneratorAction(XrayFluoDetectorConstruction* DC)
  :Detector(DC),rndmFlag("off"),randomizePrimary("off"),sigmaAngle(0),sigmaMomentum(0)
{ sigmaAngle =  1. * deg;
 sigmaMomentum = 1. * (kg*m/s);
 //default values for the spread in energy and position of the beam
 
 
 G4int n_particle = 1;
 particleGun  = new G4ParticleGun(n_particle);
 
 //create a messenger for this class
 gunMessenger = new XrayFluoPrimaryGeneratorMessenger(this);
 
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

XrayFluoPrimaryGeneratorAction::~XrayFluoPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
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

void XrayFluoPrimaryGeneratorAction::SetSigmaMomentum(G4double val)
{
  sigmaMomentum = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XrayFluoPrimaryGeneratorAction::SetSigmaAngle(G4double val)
{
  sigmaAngle = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....










