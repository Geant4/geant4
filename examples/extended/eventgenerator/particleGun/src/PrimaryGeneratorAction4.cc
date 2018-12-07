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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction4.cc
/// \brief Implementation of the PrimaryGeneratorAction4 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction4.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction4::PrimaryGeneratorAction4(G4ParticleGun* gun)
: fParticleGun(gun)
{
  // vertex volume
  //  
  G4double Rmin = 2.*mm; 
  G4double Rmax = 8.*mm;
  fRmin3 = Rmin*Rmin*Rmin;
  fRmax3 = Rmax*Rmax*Rmax;
  
  //opening angle
  //
  G4double alphaMin =  0.*deg;
  G4double alphaMax = 60.*deg;
  fCosAlphaMin = std::cos(alphaMin);
  fCosAlphaMax = std::cos(alphaMax);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction4::~PrimaryGeneratorAction4()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction4::GeneratePrimaries(G4Event* anEvent)
{  
  //vertex position uniform in spherical shell
  //
  G4double cosTheta = 2*G4UniformRand() - 1;  //cosTheta uniform in [0, pi]
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double phi      = twopi*G4UniformRand();  //phi uniform in [0, 2*pi]
  G4ThreeVector ur(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
  
  G4double R3 = fRmin3 + G4UniformRand()*(fRmax3 - fRmin3);
  G4double R  = std::pow(R3, 1./3);  
        
  fParticleGun->SetParticlePosition(R*ur);

  //particle direction uniform around ur 
  //    
  //1- in World frame
  //cosAlpha uniform in [cos(alphaMin), cos(alphaMax)]
  G4double cosAlpha = fCosAlphaMin-G4UniformRand()*(fCosAlphaMin-fCosAlphaMax);
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi      = twopi*G4UniformRand();  //psi uniform in (0,2*pi)  
  G4ThreeVector dir(sinAlpha*std::cos(psi),sinAlpha*std::sin(psi),cosAlpha);
  
  //2- rotate dir   (rotateUz transforms uz to ur)
  dir.rotateUz(ur);           

  fParticleGun->SetParticleMomentumDirection(dir);
  
  //energy
  //  
  fParticleGun->SetParticleEnergy(1*MeV);
  
  //create vertex
  //   
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
