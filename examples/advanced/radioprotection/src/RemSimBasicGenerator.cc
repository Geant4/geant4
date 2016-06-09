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
//    **********************************
//    *                                *
//    *    RemSimBasicGenerator.hh     *
//    *                                *
//    **********************************
//
// $Id: RemSimBasicGenerator.cc,v 1.5 2004/05/22 12:57:06 guatelli Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 

#include "RemSimBasicGenerator.hh"
#include "RemSimVPrimaryGeneratorFactory.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"

#ifdef G4ANALYSIS_USE  
#include "RemSimAnalysisManager.hh"
#endif

RemSimBasicGenerator::RemSimBasicGenerator()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  particleGun -> SetParticleDefinition(particleTable->FindParticle(particleName));
  particleGun -> SetParticlePosition(G4ThreeVector(0., 0., -25.*m));  
  G4ThreeVector v(0.0,0.0,1.0);
  particleGun -> SetParticleEnergy(100.*MeV);
  particleGun -> SetParticleMomentumDirection(v);
}

RemSimBasicGenerator::~RemSimBasicGenerator()
{
  delete particleGun;
}

void RemSimBasicGenerator::GeneratePrimaries(G4Event* anEvent)
{

#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun -> GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution(energy/MeV);
#endif
  particleGun -> GeneratePrimaryVertex(anEvent);
}

G4double RemSimBasicGenerator:: GetInitialEnergy()
{
 G4double primaryParticleEnergy = particleGun -> GetParticleEnergy();
 return primaryParticleEnergy;
}
