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
// $Id: HadrontherapyIonStandard.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyIonStandard.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4StepLimiter.hh"

HadrontherapyIonStandard::HadrontherapyIonStandard(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyIonStandard::~HadrontherapyIonStandard()
{ }

void HadrontherapyIonStandard::ConstructProcess()
{
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
      G4double charge = particle -> GetPDGCharge();
  
      if (particleName == "GenericIon"|| 
	  particleName == "alpha" ||
	  particleName == "He3")
	{
	  G4ionIonisation* ionisation = new G4ionIonisation();
	  G4VProcess*  multipleScattering = new G4MultipleScattering(); 
	  manager -> AddProcess(multipleScattering, -1,1,1);   
	  manager -> AddProcess(ionisation, -1,2,2);
	  manager -> AddProcess(new G4StepLimiter(),-1,-1, 3);
	}
      else
	{
	  //protons, triton, deuteron, pions and other hadrons
	  if (( charge != 0. ) && particleName != "e+" && particleName != "mu+" &&
	      particleName != "e-" && particleName != "mu-") 
	    {
	      if((!particle -> IsShortLived()) &&
		 (particle -> GetParticleName() != "chargedgeantino"))
		{
		  G4hIonisation* ionisation = new G4hIonisation();
		  G4VProcess*  multipleScattering = new G4MultipleScattering(); 
		  manager -> AddProcess(multipleScattering, -1,1,1);   
		  manager -> AddProcess(ionisation, -1,2,2);
		  manager -> AddProcess(new G4StepLimiter(),-1,-1, 3);
		}
	    }
	}
    }
}

   
