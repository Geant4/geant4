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
//
// $Id: PrimaryGeneratorAction1.cc,v 1.2 2010-07-16 07:37:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PrimaryGeneratorAction1.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction1::PrimaryGeneratorAction1(G4ParticleGun* gun)
: particleGun(gun)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction1::~PrimaryGeneratorAction1()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction1::GeneratePrimaries(G4Event* anEvent)
{
  const G4double r = 2*mm;
  const G4double zmax = 8*mm;   
  
  //vertex 1 uniform on cylinder
  //
  G4double alpha = twopi*G4UniformRand();	//alpha uniform in (0, 2*pi)
  G4double ux = std::cos(alpha);
  G4double uy = std::sin(alpha);
  G4double z = zmax*(2*G4UniformRand() - 1);	//z uniform in (-zmax, +zmax)
        
  particleGun->SetParticlePosition(G4ThreeVector(r*ux,r*uy,z));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  particleGun->SetParticleEnergy(1*MeV);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  //vertex 2 at opposite
  //
  alpha += pi;
  ux = std::cos(alpha);
  uy = std::sin(alpha);        
  particleGun->SetParticlePosition(G4ThreeVector(r*ux,r*uy,z));
  
  //particle 2 at vertex 2
  //
  const G4double dalpha = 10*deg;
  ux = std::cos(alpha + dalpha);
  uy = std::sin(alpha + dalpha);        
  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  particleGun->SetParticleEnergy(1*keV);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  //particle 3 at vertex 2
  //
  ux = std::cos(alpha - dalpha);
  uy = std::sin(alpha - dalpha);        
  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,0));    
  particleGun->SetParticleEnergy(1*GeV);
  particleGun->GeneratePrimaryVertex(anEvent);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
