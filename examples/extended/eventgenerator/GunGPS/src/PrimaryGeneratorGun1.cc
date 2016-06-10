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
/// \file eventgenerator/Gun/src/PrimaryGeneratorGun1.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorGun1.cc 68734 2013-04-05 09:47:02Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorGun1.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorGun1::PrimaryGeneratorGun1()
 : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
  // default particle kinematic
  //
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
    
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(1*eV);          
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorGun1::~PrimaryGeneratorGun1()
{
  delete fParticleGun;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorGun1::GeneratePrimaries(G4Event* anEvent)
{
    //this function is called at the begining of event
    //
    //distribution uniform in solid angle
    //
    G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
    G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    G4double ux = sinTheta*std::cos(phi),
             uy = sinTheta*std::sin(phi),
             uz = cosTheta;

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
