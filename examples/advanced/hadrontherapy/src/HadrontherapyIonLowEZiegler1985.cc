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
// $Id: HadrontherapyIonLowEZiegler1985.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, G. Candiano, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// --------------------------------------------------------------

#include "HadrontherapyIonLowEZiegler1985.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4StepLimiter.hh"
#include "G4hZiegler1985p.hh"

HadrontherapyIonLowEZiegler1985::HadrontherapyIonLowEZiegler1985(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyIonLowEZiegler1985::~HadrontherapyIonLowEZiegler1985()
{ }

void HadrontherapyIonLowEZiegler1985::ConstructProcess()
{
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
      G4double charge = particle -> GetPDGCharge();

      // protons
      if (particleName == "proton")
	{
         G4hLowEnergyIonisation* ionisation = new G4hLowEnergyIonisation();   
         // Electronic Stopping Power: Ziegler 1985 parameterisation
         // Nuclear stopping power: Ziegler 1985 parameterisation   
         ionisation -> SetElectronicStoppingPowerModel(particle, "Ziegler1985p");
         ionisation -> SetNuclearStoppingPowerModel("Ziegler1985");
	 ionisation -> SetNuclearStoppingOn() ; 

         G4VProcess*  multipleScattering = new G4MultipleScattering(); 
         manager -> AddProcess(multipleScattering, -1,1,1);   
         manager -> AddProcess(ionisation, -1,2,2);
	 manager -> AddProcess(new G4StepLimiter(),-1,-1, 3);
	}

      // ions, alpha, pions, kaons, generic ions.....
    if (( charge != 0. ) && (particleName != "e+") && (particleName != "mu+") &&
	  (particleName != "e-") && (particleName != "mu-") && (particleName != "proton")) 
	{
	  if((!particle -> IsShortLived()) &&
	     (particle -> GetParticleName() != "chargedgeantino"))
	    {
	      G4hLowEnergyIonisation* ionisation = new G4hLowEnergyIonisation();
		
	      ionisation -> SetNuclearStoppingOn() ;
		 
	      G4VProcess*  multipleScattering = new G4MultipleScattering(); 
	      manager -> AddProcess(multipleScattering, -1,1,1);   
	      manager -> AddProcess(ionisation, -1,2,2);
	      manager -> AddProcess(new G4StepLimiter(),-1,-1, 3);
	    }
	}
    }
}

   
