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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2009/03/21 18:37:27 vnivanch Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"
#include "G4Geantino.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* theDete) : 
  theDetector(theDete)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle); 

  G4ParticleDefinition* particle
    = G4Gamma::GammaDefinition();
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(1.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(-5.*cm,0.*cm,0*cm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::GetPrimaryEnergy()
{
  return particleGun->GetParticleEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (particleGun->GetParticleDefinition() == G4Geantino::GeantinoDefinition())
    {
      particleGun->GeneratePrimaryVertex(anEvent);
      return;
    }

  G4int detType = theDetector->GetGeometryType();
  
  //Isotropic direction
  G4double cosTheta = -1.0 + 2.0*G4UniformRand();
  G4double sinTheta = std::sqrt(1.0-cosTheta*cosTheta);
  G4double phi = twopi*G4UniformRand();
  G4ThreeVector direction(sinTheta*std::sin(phi),sinTheta*std::cos(phi),
			  cosTheta);
  particleGun->SetParticleMomentumDirection(direction);

  if (detType == 2) //point-like source, 5 mm above the casing
    {
      G4double casingHeight = theDetector->GetCasingHeight();
      G4double sourceToCasing = 5.0*mm;
      G4double zpos = casingHeight/2. + sourceToCasing;
      particleGun->SetParticlePosition(G4ThreeVector(0,0,zpos));

    }
  else if (detType == 3) //uniformly sample from the water block
    {
      G4double WaterZpos = theDetector->GetWaterZPosition();
      G4double waterHeight = theDetector->GetWaterHeight();
      G4double waterRadius = theDetector->GetWaterDiameter()/2.;
   
      G4double b     = waterRadius*waterRadius;
      G4double r     = std::sqrt(b*G4UniformRand());
      G4double phi   = twopi * G4UniformRand();

      G4double x = r*std::cos(phi);
      G4double y = r*std::sin(phi);
      G4double z = (G4UniformRand()-0.5)*waterHeight;

      G4ThreeVector initialPos(x,y,z+WaterZpos);
      particleGun->SetParticlePosition(initialPos);
    }

  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

