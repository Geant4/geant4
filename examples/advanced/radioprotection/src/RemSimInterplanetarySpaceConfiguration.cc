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

#include "CLHEP/Random/RandGeneral.h"
#include "RemSimInterplanetarySpaceConfiguration.hh"
#include "RemSimVPrimaryGeneratorFactory.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#ifdef G4ANALYSIS_USE  
#include "RemSimAnalysisManager.hh"
#endif
#include "RemSimRunAction.hh"
#include <fstream>
#include <strstream>

RemSimInterplanetarySpaceConfiguration::RemSimInterplanetarySpaceConfiguration()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  run = new RemSimRunAction();
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  particleGun -> SetParticleDefinition(particleTable->FindParticle(particleName));
  particleGun -> SetParticlePosition(G4ThreeVector(0., 0., -25.*m));  
  G4ThreeVector v(0.0,0.0,1.0);
  particleGun -> SetParticleEnergy(1.*MeV);
  particleGun -> SetParticleMomentumDirection(v);
}

RemSimInterplanetarySpaceConfiguration::~RemSimInterplanetarySpaceConfiguration()
{
  delete run;
  delete particleGun;
}

void RemSimInterplanetarySpaceConfiguration::GeneratePrimaries(G4Event* anEvent){
  // Read the ASCII files containing  the fluxes of particles in respect
  // to the energy in MeV

  G4DataVector* energies = run -> GetPrimaryParticleEnergy();
  G4DataVector* data = run -> GetPrimaryParticleEnergyDistribution();	 
  G4double sum = run -> GetPrimaryParticleEnergyDistributionSum();
  G4double partSum = 0;
  G4int j = 0;
  G4double random = sum*G4UniformRand();
  while (partSum<random)
    {
      partSum += (*data)[j];
      j++;
    }
	 
  particleGun -> SetParticleEnergy((*energies)[j]);    
 
#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun -> GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution(energy/MeV);
#endif
 
 particleGun -> GeneratePrimaryVertex(anEvent);
}

G4double RemSimInterplanetarySpaceConfiguration:: GetInitialEnergy()
{
 G4double primaryParticleEnergy = particleGun -> GetParticleEnergy();
 return primaryParticleEnergy;
}
