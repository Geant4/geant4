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
// $Id: RemSimPrimaryGeneratorAction.cc,v 1.17 2006-06-29 16:24:09 gunter Exp $// Author: Susanna Guatelli, guatelli@ge.infn.it

#include <fstream>
#include <cmath>

#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimPrimaryGeneratorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"

RemSimPrimaryGeneratorAction::RemSimPrimaryGeneratorAction()
{
  readFile = false;
 
  particleGun = new G4ParticleGun(1);
  messenger = new RemSimPrimaryGeneratorMessenger(this);
  
  energies = new G4DataVector;
  data = new G4DataVector;
 
  size = 0;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String ParticleName = "proton";
  G4ParticleDefinition* particle = particleTable->FindParticle(ParticleName);

  particleGun->SetParticleDefinition(particle); 

  // Default values of the primary generator
  G4ThreeVector position(0,0,-25. *m); 
  particleGun -> SetParticlePosition(position);
  G4ThreeVector direction(0., 0., 1.);
  particleGun -> SetParticleMomentumDirection(direction);
  particleGun -> SetParticleEnergy(250. * MeV); 
}

RemSimPrimaryGeneratorAction::~RemSimPrimaryGeneratorAction()
{
  delete data;
  delete energies;
  delete [] cumulate; 
  delete [] energy;
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
 G4int baryon = particleGun -> GetParticleDefinition() -> GetBaryonNumber();

 // Generate the energy spectrum of primary particles according to 
 // the flux reported in the .txt file (firt column = Energy (MeV/nucl),
 // second column = differnetila flux). The fluxes are derived from CREME96  
 if (readFile == true)
  { 
  // Uniform number between 0. and 1.
  G4double random = G4UniformRand();
  
  G4int nbelow = 0;	  // largest k such that I[k] is known to be <= rand
  G4int nabove = size;  // largest k such that I[k] is known to be >  rand
  G4int middle;
  
  while (nabove > nbelow+1) {
    middle = (nabove + nbelow+1)>>1;
    if (random >= cumulate[middle]) {
      nbelow = middle;
    } else {
      nabove = middle;
    }
  }
   
 // Linear interpolation

  G4double energy_gun = 0. * MeV;
 
  G4double binMeasure = cumulate[nabove] - cumulate[nbelow];

    if ( binMeasure == 0 ) { 

      energy_gun = (energy[nbelow] + random*(energy[nabove] - energy[nbelow]))* baryon;
      G4cout<<"BinMeasure is zero !!!!!"<<G4endl;
    }

    else
      {
    G4double binFraction = (random - cumulate[nbelow]) / binMeasure;
    energy_gun = energy[nbelow] + binFraction *(energy[nabove] - energy[nbelow]);
      } 
  particleGun -> SetParticleEnergy(energy_gun * baryon); 
  }


particleGun -> GeneratePrimaryVertex(anEvent);
}

void RemSimPrimaryGeneratorAction::ReadProbability(G4String fileName)
{
  // This method allows to read the .txt files containing the
  // fluxes of the galactic cosmic rays derived form CREME96
  readFile = true;
  
  std::ifstream file(fileName, std::ios::in);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
	  G4String excep = "RemSimPrimaryGenerator - data file: " 
                            + fileName + " not found";
	  G4Exception("RemSimPrimaryGeneratorAction::ReadProbability()",
		      "RadioP001",FatalException,excep);
    }

  G4double a = 0;
  G4int index = 0;
  G4int k = 1;
 
 do
    {
      file >> a; 
      G4int nColumns = 2;
      // The file is organized into two columns:
      // column contains the probability distribution
      // The file terminates with the pattern: -1   
      //                                       -2   
      if (a == -1 || a == -2)
	{
        
	}
      else
	{
          if (k%nColumns != 0)
	    {	
	      G4double ene = a *MeV;
	      energies -> push_back(ene);
	      //   G4cout << "energy: " << ene << G4endl;
              k++;
	    }
	  else if (k%nColumns == 0)
	    {
	   G4double data_value = a;
           data -> push_back(data_value);
          //G4cout<< "probability: "<< data_value << G4endl; 
	   k=1;
           index++; 
	    }  
	}	 
    } while (a != -2); // end of file
  
  file.close();
  size = index -1 ;
  cumulate = new G4double[size+1];
  energy = new G4double[size+1];
  cumulate[0] = 0;
  energy [0] = (*energies)[0];

  for (G4int ptn = 0; ptn <size; ptn++ ) {
 
  G4double yy = ((*data)[ptn+1] + (*data)[ptn])/2.;
  G4double weight = yy *((*energies)[ptn +1] - (*energies)[ptn]);
  cumulate[ptn +1] = cumulate[ptn] + weight;
  energy[ptn+1] = (*energies)[ptn+1];
 } 
  
 // Normalise the cumulate to 1 (that corresponds to the last value)
 for ( G4int ptn= 0; ptn < size + 1; ptn++ ) 
{            
  cumulate[ptn] /= cumulate[size];  
}
}


