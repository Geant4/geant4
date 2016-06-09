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
// ********************************************************************

//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorActionIr.cc     *
//    *                                          *
//    ********************************************
//
// $Id: BrachyPrimaryGeneratorActionIr.cc,v 1.6 2003/05/22 17:20:43 guatelli Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
#include "BrachyPrimaryGeneratorActionIr.hh"

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

BrachyPrimaryGeneratorActionIr::BrachyPrimaryGeneratorActionIr()
{
  G4int NumParticles = 1;
  particleGun = new G4ParticleGun(NumParticles); 
}

BrachyPrimaryGeneratorActionIr::~BrachyPrimaryGeneratorActionIr()
{
  if(particleGun)
    delete particleGun;
}

void BrachyPrimaryGeneratorActionIr::GeneratePrimaries(G4Event* anEvent)
{
#ifdef G4ANALYSIS_USE
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
#endif

  G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "gamma";
  G4ParticleDefinition* pParticle = pParticleTable->FindParticle(ParticleName);
  particleGun->SetParticleDefinition(pParticle);
 
  //  Random generation of gamma source point inside the Iodium core
  G4double x,y,z;
  G4double radius= 0.30*mm;
  do{
    x = (G4UniformRand()-0.5)*(radius)/0.5;
    y = (G4UniformRand()-0.5)*(radius)/0.5;
  }while(x*x+y*y > radius*radius);
 
  z = (G4UniformRand()-0.5)*1.75*mm/0.5 ;

  G4ThreeVector position(x,y,z);
  particleGun->SetParticlePosition(position);
  // Random generation of the impulse direction
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

  primaryParticleEnergy = 356*keV;
  particleGun->SetParticleEnergy(primaryParticleEnergy);
 
  //1D Histogram of primary particle energy ...
#ifdef G4ANALYSIS_USE
  analysis->PrimaryParticleEnergySpectrum(primaryParticleEnergy);
#endif   
  particleGun->GeneratePrimaryVertex(anEvent);
}

