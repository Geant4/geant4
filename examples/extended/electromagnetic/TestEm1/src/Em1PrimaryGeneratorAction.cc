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
//
// $Id: Em1PrimaryGeneratorAction.cc,v 1.5 2001-12-07 11:49:09 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em1PrimaryGeneratorAction.hh"

#include "Em1DetectorConstruction.hh"
#include "Em1PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PrimaryGeneratorAction::Em1PrimaryGeneratorAction(
                                               Em1DetectorConstruction* Em1DC)
:Em1Detector(Em1DC)					       
{
  particleGun  = new G4ParticleGun(1);
  SetDefaultKinematic(0);
  rndmBeam = 0.;
    
  //create a messenger for this class
  gunMessenger = new Em1PrimaryGeneratorMessenger(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em1PrimaryGeneratorAction::~Em1PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PrimaryGeneratorAction::SetDefaultKinematic(G4int front)
{
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(100*MeV);
  G4double position = 0.*cm;
  if (front) position = -0.5*(Em1Detector->GetSize());
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  //randomize the beam, if requested.
  if (rndmBeam > 0.) 
    {
      G4ThreeVector oldPosition = particleGun->GetParticlePosition();    
      G4double rbeam = 0.5*(Em1Detector->GetSize())*rndmBeam;
      G4double x0 = oldPosition.x();
      G4double y0 = oldPosition.y() + (2*G4UniformRand()-1.)*rbeam;
      G4double z0 = oldPosition.z() + (2*G4UniformRand()-1.)*rbeam;
      particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticlePosition(oldPosition);      
    }
  else  particleGun->GeneratePrimaryVertex(anEvent); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

