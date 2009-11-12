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
// $Id: RemSimIonICRU.cc,v 1.11 2009-11-12 05:12:18 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimIonICRU.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4StepLimiter.hh"
#include "G4IonParametrisedLossModel.hh"

RemSimIonICRU::RemSimIonICRU(const G4String& name): G4VPhysicsConstructor(name)
{ }

RemSimIonICRU::~RemSimIonICRU()
{ }

void RemSimIonICRU::ConstructProcess()
{
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
      G4double charge = particle -> GetPDGCharge();

      if (particleName  == "proton" 
	  || particleName == "pi+"
          || particleName == "pi-")
	{
	  manager -> AddProcess(new G4hMultipleScattering, -1,1,1);   

	  G4hIonisation* ionisation = new G4hIonisation();
	  ionisation -> SetStepFunction(0.2, 50*um);
	  manager -> AddProcess(ionisation, -1,2,2); 
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}

      else if(particleName == "alpha"      ||
	      particleName == "deuteron"   ||
	      particleName == "triton"     ||
	      particleName == "He3")
	{
	  //multiple scattering
	  manager->AddProcess(new G4hMultipleScattering,-1,1,1);
	
	  //ionisation
	  G4ionIonisation* ionIoni = new G4ionIonisation();
	  ionIoni->SetStepFunction(0.1, 20*um);
	  manager->AddProcess(ionIoni,-1, 2, 2);
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}

      else if (particleName == "GenericIon")
	{
	  // OBJECT may be dynamically created as either a GenericIon or nucleus
	  // G4Nucleus exists and therefore has particle type nucleus
	  // genericIon:
	
	  //multiple scattering
	  manager->AddProcess(new G4hMultipleScattering,-1,1,1);

	  //ionisation
	  G4ionIonisation* ionIoni = new G4ionIonisation();
	  ionIoni->SetEmModel(new G4IonParametrisedLossModel());
	  ionIoni->SetStepFunction(0.1, 20*um);
	  manager->AddProcess(ionIoni,-1, 2, 2);
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}
      else{     
	if (( charge != 0. ) && particleName != "e+" && particleName != "mu+" &&
            particleName != "e-" && particleName != "mu-")
	  {
	    if((!particle -> IsShortLived()) &&
	       (particle -> GetParticleName() != "chargedgeantino"))
	      {
		G4hIonisation* ionisation = new G4hIonisation();
		G4VProcess*  multipleScattering = new G4hMultipleScattering(); 
		manager -> AddProcess(multipleScattering, -1,1,1);   
		manager -> AddProcess(ionisation, -1,2,2);
		manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	      }
	  }
      }
    }
}

   
