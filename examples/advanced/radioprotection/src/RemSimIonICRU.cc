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
// $Id: RemSimIonICRU.cc,v 1.10 2006/06/29 16:23:51 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimIonICRU.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4StepLimiter.hh"

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

      if (particleName  == "proton" || particleName == "GenericIon"
          || particleName == "alpha" )
	{
	  G4hLowEnergyIonisation* ionisation = new G4hLowEnergyIonisation();
	  G4VProcess*  multipleScattering = new G4MultipleScattering(); 
	  manager -> AddProcess(multipleScattering, -1,1,1);   
	  manager -> AddProcess(ionisation, -1,2,2);
	  manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}
      else{     
	if (( charge != 0. ) && particleName != "e+" && particleName != "mu+" &&
            particleName != "e-" && particleName != "mu-") 
	  {
	    if((!particle -> IsShortLived()) &&
	       (particle -> GetParticleName() != "chargedgeantino"))
	      {
		G4hLowEnergyIonisation* ionisation = new G4hLowEnergyIonisation();
		G4VProcess*  multipleScattering = new G4MultipleScattering(); 
		manager -> AddProcess(multipleScattering, -1,1,1);   
		manager -> AddProcess(ionisation, -1,2,2);
		manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	      }
	  }
      }}
}

   
