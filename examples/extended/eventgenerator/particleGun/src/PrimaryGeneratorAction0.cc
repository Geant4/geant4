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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction0.cc
/// \brief Implementation of the PrimaryGeneratorAction0 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction0.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction0::PrimaryGeneratorAction0(G4ParticleGun* gun)
: fParticleGun(gun)
{
  //solid angle
  //
  G4double alphaMin =  0*deg;      //alpha in [0,pi]
  G4double alphaMax = 60*deg;
  fCosAlphaMin = std::cos(alphaMin);
  fCosAlphaMax = std::cos(alphaMax);
  
  fPsiMin = 0*deg;       //psi in [0, 2*pi]
  fPsiMax = 360*deg;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction0::~PrimaryGeneratorAction0()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction0::GeneratePrimaries(G4Event* anEvent)
{  
  //vertex position fixed
  //
  G4double x0 = 0*mm, y0 = 0*mm, z0 = 0*mm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  //direction uniform in solid angle
  //
  G4double cosAlpha = fCosAlphaMin-G4UniformRand()*(fCosAlphaMin-fCosAlphaMax);
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi = fPsiMin + G4UniformRand()*(fPsiMax - fPsiMin);

  G4double ux = sinAlpha*std::cos(psi),
           uy = sinAlpha*std::sin(psi),
           uz = cosAlpha;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
