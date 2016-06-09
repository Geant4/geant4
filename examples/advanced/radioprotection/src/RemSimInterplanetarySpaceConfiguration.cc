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
// Author : Susanna Guatelli, guatelli@ge.infn.it
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

  moon = false;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  particleGun -> SetParticleDefinition(particleTable->FindParticle(particleName));
  particleGun -> SetParticlePosition(G4ThreeVector(0., 0., -25.*m));  
  G4ThreeVector v(0.0,0.0,1.0);
  particleGun -> SetParticleEnergy(1.*MeV);
  particleGun -> SetParticleMomentumDirection(v);
 
  energies = new G4DataVector;
  data = new G4DataVector;
}

RemSimInterplanetarySpaceConfiguration::~RemSimInterplanetarySpaceConfiguration()
{
  delete data;
  delete energies;
  delete particleGun;
}

void RemSimInterplanetarySpaceConfiguration::GeneratePrimaries(G4Event* anEvent){
  G4double sum = GetPrimaryParticleEnergyDistributionSum();
  G4double partSum = 0;
  G4int j = 0;
  G4double random = sum*G4UniformRand();
  while (partSum<random)
    {
      partSum += (*data)[j];
      j++;
    }

    G4String particleName = particleGun -> GetParticleDefinition() -> 
  GetParticleName();	 
  
  G4int n = 0;

  if (particleName == "alpha") n = 4;
  else {
    if (particleName == "IonC12") n = 12;
    else
      { if (particleName == "IonSi28") n = 28;
      else {
        if (particleName == "IonFe52") n = 52;
	else {
	  if (particleName == "IonO16") n = 16;
	  else {
	    if (particleName == "proton") n = 1 ;
	  }
	}
      }
      }
  }

  particleGun -> SetParticleEnergy((*energies)[j] * n);    
 
#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun -> GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution((energy/n) / MeV);
#endif
 
 if (moon == true) MoonConfiguration();
 else
   {
     particleGun -> SetParticlePosition(G4ThreeVector(0., 0., -25.*m));  
     G4ThreeVector v(0.0,0.0,1.0); 
     particleGun -> SetParticleMomentumDirection(v);
   }
  
 particleGun -> GeneratePrimaryVertex(anEvent);
}

void RemSimInterplanetarySpaceConfiguration:: MoonConfiguration() 
{
 //Generate the primary particles on a hemisphere with random direction
 //position
  G4double radius = 25.* m;
  G4double angle = pi * G4UniformRand()*rad;
  G4double y0 = radius*std::cos(angle);
  G4double x0 = 0.*m;
  G4double z0 = -radius*std::sin(angle);

  if ( z0 < 0. *m)
    {
   particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    }

 //direction 
      G4double angledir = pi * G4UniformRand()*rad;
      if ((angledir> (pi/2.-angle)) && angledir<(3*pi/2.-angle))
	{
      G4double a,b,c;
      b=std::cos(angledir);
      a=0.;
      c=std::sin(angledir);
      G4ThreeVector direction(a,b,c);
      particleGun -> SetParticleMomentumDirection(direction);
	}
}

G4double RemSimInterplanetarySpaceConfiguration:: GetInitialEnergy()
{
 G4double primaryParticleEnergy = particleGun -> GetParticleEnergy();
 return primaryParticleEnergy;
}
void RemSimInterplanetarySpaceConfiguration::SetMoon(G4bool value)
{
  moon = value;
}

void RemSimInterplanetarySpaceConfiguration::Read(G4String name)
{    
  ReadData(MeV,name);
  G4cout << name << "  is the input file!" << G4endl;
}

void RemSimInterplanetarySpaceConfiguration::ReadData(G4double unitE, G4String fileName)
{
  char nameChar[100] = {""};
  std::ostrstream ost(nameChar, 100, std::ios::out);
 
  ost << fileName;
  
  G4String name(nameChar);
  
  std::ifstream file(fileName);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
	  G4String excep = "RemSimInterplanteraySpaceConfiguration - data file: not found";
	  G4Exception(excep);
    }
  G4double a = 0;
  G4int k = 1;
  
  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	  
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unitE;
	      energies->push_back(e);  
	      //              G4cout<<e<<"energy";
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a;
	      data->push_back(value);
	      //G4cout<<" "<<a<<"flux"<<G4endl;
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
}
G4double RemSimInterplanetarySpaceConfiguration::GetPrimaryParticleEnergyDistributionSum()
{
  G4double sum = 0;
  size_t size = data -> size();
  for (size_t i = 0; i < size; i++)
     {
       sum+=(*data)[i];
     }
   return sum;
 }
