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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.12 2005-09-08 06:56:18 guatelli Exp $// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimPrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "globals.hh"
#ifdef G4ANALYSIS_USE  
#include "RemSimAnalysisManager.hh"
#endif
#include <fstream>
#include <strstream>
#include "Randomize.hh"

RemSimPrimaryGeneratorAction::RemSimPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  //moon = false;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName = "proton";
  
  particleGun -> SetParticleDefinition(particleTable->FindParticle(particleName));
  particleGun -> SetParticlePosition(G4ThreeVector(0., 0., -25.*m));  
  G4ThreeVector v(0.0,0.0,1.0);
  particleGun -> SetParticleEnergy(1.*GeV);
  particleGun -> SetParticleMomentumDirection(v);
  value = "Basic";
  messenger = new RemSimPrimaryGeneratorMessenger(this);
  energies = new G4DataVector;
  data = new G4DataVector;
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
  delete data;
  delete energies;
  delete messenger;
  delete particleGun;
}

G4double RemSimPrimaryGeneratorAction::GetInitialEnergy()
{
   G4double primaryParticleEnergy = particleGun -> GetParticleEnergy();
   return primaryParticleEnergy;   
}

void RemSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4int n = 1;

  if(value == "Interplanetary") 
   { 
     G4double sum = GetPrimaryParticleEnergyDistributionSum();
     G4double partSum = 0;
     G4int j = 0;
     G4double random = sum*G4UniformRand();
     while (partSum<random)
       {
	 partSum += (*data)[j];
	 j++;
       }

     G4String particleName = particleGun -> GetParticleDefinition() -> GetParticleName();	 
     if (particleName == "alpha") n = 4;
     else{ if (particleName == "proton") n = 1 ;}
     particleGun -> SetParticleEnergy((*energies)[j] * n);    
 }
 
#ifdef G4ANALYSIS_USE   
 G4double energy = particleGun -> GetParticleEnergy(); 
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> primaryParticleEnergyDistribution((energy/n) / MeV);
#endif

particleGun -> GeneratePrimaryVertex(anEvent);
}

void RemSimPrimaryGeneratorAction::SelectPrimaries(G4String val)
{ 
  value = val;

  if(value == "Basic") 
    G4cout<< "The configuration is the basic generator" <<G4endl;

  else if(value == "Interplanetary") 
    {
      G4cout<< "The configuration is the interplanetary space configuration" 	    <<G4endl;
     
      G4cout<< 
	"Remember to type /run/data file.txt with the energy spectrum of primary particles!!!" <<G4endl;
     
    }
  else G4cout << "This Generator is not defined!" <<G4endl;  
}

void RemSimPrimaryGeneratorAction::Read(G4String fileName)
{
  ReadData(MeV,fileName);
  G4cout << fileName << "  is the input file!" << G4endl;
}

void RemSimPrimaryGeneratorAction::ReadData(G4double unitE, G4String fileName)
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
	      
	      k++;
	      
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a;
	      data->push_back(value);
	      
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
}
G4double RemSimPrimaryGeneratorAction::GetPrimaryParticleEnergyDistributionSum()
{
  G4double sum = 0;
  size_t size = data -> size();
  for (size_t i = 0; i < size; i++)
     {
       sum+=(*data)[i];
     }
   return sum;
 }
