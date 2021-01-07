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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction3.cc
/// \brief Implementation of the PrimaryGeneratorAction3 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PrimaryGeneratorAction3.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction3::PrimaryGeneratorAction3(G4ParticleGun* gun)
: fParticleGun(gun)
{    
  // direction
  //
  G4double theta = 90*deg, phi = 45*deg;
  fNewUz = G4ThreeVector(std::sin(theta)*std::cos(phi),
                         std::sin(theta)*std::sin(phi),
                         std::cos(theta));
                        
  fAlphaMax = 15*deg;                        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction3::~PrimaryGeneratorAction3()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction3::GeneratePrimaries(G4Event* anEvent)
{
  //particle direction uniform around fNewUz axis
  //
  //1- in World frame      
  //cosAlpha uniform in [cos(0), cos(fAlphaMax)]
  G4double cosAlpha = 1. - G4UniformRand()*(1.- std::cos(fAlphaMax));
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi      = twopi*G4UniformRand();  //psi uniform in [0, 2*pi]  
  G4ThreeVector dir(sinAlpha*std::cos(psi),sinAlpha*std::sin(psi),cosAlpha);
  
  //2- rotate dir   (rotateUz transforms uz to fNewUz)
  dir.rotateUz(fNewUz);           

  fParticleGun->SetParticleMomentumDirection(dir);
  
  //set energy
  //
  fParticleGun->SetParticleEnergy(1*MeV);    

  //create vertex
  //   
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
