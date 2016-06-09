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
// $Id: BrachyPrimaryGeneratorActionI.cc,v 1.7 2003/12/09 15:30:01 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
#include "BrachyPrimaryGeneratorActionI.hh"

#ifdef G4ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "globals.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

BrachyPrimaryGeneratorActionI::BrachyPrimaryGeneratorActionI()
{
  G4int numberParticles = 1;
 
  particleGun = new G4ParticleGun(numberParticles);
  
  // Gamma energy spectrum ...
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
#ifdef G4ANALYSIS_USE
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
#endif

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "gamma";
  G4ParticleDefinition* particle = particleTable->FindParticle(ParticleName);

  particleGun->SetParticleDefinition(particle);

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
  particleGun->SetParticlePosition(position);

  // Random generation of the impulse direction of primary particles ...
  G4double a,b,c;
  G4double n;
  do{
    a = (G4UniformRand()-0.5)/0.5;
    b = (G4UniformRand()-0.5)/0.5; 
    c = (G4UniformRand()-0.5)/0.5;
    n = a*a+b*b+c*c;
  }while(n > 1 || n == 0.0);
  n = sqrt(n);
  a /= n;
  b /= n;
  c /= n;

  G4ThreeVector direction(a,b,c);
  particleGun->SetParticleMomentumDirection(direction);
 
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

  particleGun->SetParticleEnergy(primaryParticleEnergy);

  // Fill 1D histogram with gamma energy spectrum ...
#ifdef G4ANALYSIS_USE    
  analysis-> PrimaryParticleEnergySpectrum(primaryParticleEnergy);
#endif  
  // generate primary particle
  particleGun->GeneratePrimaryVertex(anEvent);
}
 
G4double  BrachyPrimaryGeneratorActionI::GetEnergy()
{ 
  return primaryParticleEnergy;
}

