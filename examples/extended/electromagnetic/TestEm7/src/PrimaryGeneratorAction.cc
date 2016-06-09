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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2003/04/22 16:25:07 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:detector(det)					       
{
  particleGun  = new G4ParticleGun(1);
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleEnergy(160*MeV);  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
    
  rndmBeam   = 0.;
  EbeamCumul = 0.;
    
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  G4double x0 = -0.5*(detector->GetWorldSizeX());
  G4double y0 = 0.*cm, z0 = 0.*cm;
    
  //randomize the beam, if requested.
  //
  if (rndmBeam > 0.) 
    {
      if (rndmBeam > detector->GetAbsorSizeYZ())
        rndmBeam = detector->GetAbsorSizeYZ(); 
      G4double rbeam = 0.5*rndmBeam;
      y0 = (2*G4UniformRand()-1.)*rbeam;
      z0 = (2*G4UniformRand()-1.)*rbeam;
    }
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));  
  particleGun->GeneratePrimaryVertex(anEvent);
  
  EbeamCumul += particleGun->GetParticleEnergy(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

