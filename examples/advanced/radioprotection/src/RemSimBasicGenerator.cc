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
// $Id: RemSimBasicGenerator.cc,v 1.3 2004-03-12 10:55:55 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimBasicGenerator.hh"
#include "RemSimVPrimaryGeneratorFactory.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#ifdef G4ANALYSIS_USE  
#include "RemSimAnalysisManager.hh"
#endif
RemSimBasicGenerator::RemSimBasicGenerator():randomDirection("off")
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
  particleGun->SetParticlePosition(G4ThreeVector(0., 0., -5.*m));  
  G4ThreeVector v(0.0,0.0,1.0);
  particleGun->SetParticleEnergy(1.*MeV);
  particleGun->SetParticleMomentumDirection(v);
}

RemSimBasicGenerator::~RemSimBasicGenerator()
{
  delete particleGun;
}

void RemSimBasicGenerator::GeneratePrimaries(G4Event* anEvent)
{
if(randomDirection == "on")
    {
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
      particleGun -> SetParticleMomentumDirection(direction);
    }
 

#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun->GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution(energy/MeV);
#endif
  particleGun->GeneratePrimaryVertex(anEvent);


}

void RemSimBasicGenerator:: GenerateIsotropicFlux()
{
  randomDirection = "on";
}

G4double RemSimBasicGenerator:: GetInitialEnergy()
{
 G4double primaryParticleEnergy = particleGun->GetParticleEnergy();
 return primaryParticleEnergy;
}
