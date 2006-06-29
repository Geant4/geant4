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
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************************
//    *                                          *
//    *    BrachyPrimaryGeneratorActionIr.cc     *
//    *                                          *
//    ********************************************
//
// $Id: BrachyPrimaryGeneratorActionIr.cc,v 1.11 2006-06-29 15:48:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//#include "G4IonTable.hh"
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

  // Define primary particle type
  G4ParticleTable* pParticleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "gamma";
  G4ParticleDefinition* pParticle = pParticleTable -> FindParticle(ParticleName);
  particleGun -> SetParticleDefinition(pParticle);
 
  //  Random generation of gamma source point inside the Iridium core
  G4double x,y,z;
  G4double radius = 0.30*mm;
  do{
    x = (G4UniformRand()-0.5)*(radius)/0.5;
    y = (G4UniformRand()-0.5)*(radius)/0.5;
  }while(x*x+y*y > radius*radius);
 
  z = (G4UniformRand()-0.5)*1.75*mm/0.5 -1.975*mm  ;

  G4ThreeVector position(x,y,z);
  particleGun -> SetParticlePosition(position);

  // Random generation of the impulse direction
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
  particleGun->SetParticleMomentumDirection(direction);

  // Primary particle energy
  primaryParticleEnergy = 356.*keV;
  particleGun->SetParticleEnergy(primaryParticleEnergy);
 
  //1D Histogram of primary particle energy ...
#ifdef G4ANALYSIS_USE
  analysis -> PrimaryParticleEnergySpectrum(primaryParticleEnergy);
#endif   
  
  // Generate a primary particle
  particleGun -> GeneratePrimaryVertex(anEvent);
}

