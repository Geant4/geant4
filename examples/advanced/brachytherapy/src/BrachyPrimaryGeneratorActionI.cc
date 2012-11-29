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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S.Guatelli
//
//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorActionI.cc      *
//    *                                          *
//    ********************************************
//
// $Id$
//
#include "BrachyPrimaryGeneratorActionI.hh"
#include "globals.hh"
#include "Randomize.hh"  
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

BrachyPrimaryGeneratorActionI::BrachyPrimaryGeneratorActionI()
{
  G4int numberParticles = 1;
 
  particleGun = new G4ParticleGun(numberParticles);
  
  // Gamma energy spectrum ...
  // Fill a vector with the energy probabilities
  energySpectrum.push_back(0.783913);
  energySpectrum.push_back(0.170416);
  energySpectrum.push_back(0.045671); 
}

BrachyPrimaryGeneratorActionI::~BrachyPrimaryGeneratorActionI()
{
 if(particleGun)
	delete particleGun;
}

void BrachyPrimaryGeneratorActionI::GeneratePrimaries(G4Event* anEvent)
{
  // Define the primary particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "gamma";
  G4ParticleDefinition* particle = particleTable -> FindParticle(ParticleName);

  particleGun -> SetParticleDefinition(particle);

  //  Random generation of gamma source point inside the Iodium core ...
  G4double x,y,z;
  G4double radiuMax = 0.30*mm;
  G4double radiusMin = 0.085*mm;
 
  do{
     x = (G4UniformRand()-0.5)*(radiuMax)/0.5;
     y = (G4UniformRand()-0.5)*(radiuMax)/0.5;
  }while(((x*x+y*y )> (radiuMax*radiuMax))||((x*x+y*y)<(radiusMin*radiusMin)));
 
  z = (G4UniformRand()-0.5)*1.75*mm/0.5 ;

  G4ThreeVector position(x,y,z);
  particleGun -> SetParticlePosition(position);

  // Random generation of the impulse direction of primary particles ...
  G4double a,b,c;
  G4double n;
  do{
    a = (G4UniformRand()-0.5)/0.5;
    b = (G4UniformRand()-0.5)/0.5; 
    c = (G4UniformRand()-0.5)/0.5;
    n = a*a+b*b+c*c;
  }while(n > 1 || n == 0.0);
  n = std::sqrt(n);
  a /= n;
  b /= n;
  c /= n;

  G4ThreeVector direction(a,b,c);
  particleGun -> SetParticleMomentumDirection(direction);
 
  // Generate the primary particles with a defined energy spectrum
  G4double random = G4UniformRand();
  G4double sum = 0;
  G4int i = 0;
  while(sum<random){sum+=energySpectrum[i];
  i++;}

  // energy spectrum
  if(i==1){primaryParticleEnergy = 27.4*keV;}
  else{ 
    if(i==2){primaryParticleEnergy = 31.4*keV;}
    else {primaryParticleEnergy = 35.5*keV;}}

  particleGun -> SetParticleEnergy(primaryParticleEnergy);

  // generate primary particle
  particleGun->GeneratePrimaryVertex(anEvent);
}
 
G4double  BrachyPrimaryGeneratorActionI::GetEnergy()
{ 
  return primaryParticleEnergy;
}

