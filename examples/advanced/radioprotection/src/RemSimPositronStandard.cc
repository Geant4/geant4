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
// $Id: RemSimPositronStandard.cc,v 1.5 2005/09/08 06:56:18 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it
#include "RemSimPositronStandard.hh"
#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4StepLimiter.hh"

RemSimPositronStandard::RemSimPositronStandard(const G4String& name): G4VPhysicsConstructor(name)
{ }

RemSimPositronStandard::~RemSimPositronStandard()
{ }

void RemSimPositronStandard::ConstructProcess()
{
  // Add standard processes for positrons
  
  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if (particleName == "e+") 
	{
	  manager -> AddProcess(new G4MultipleScattering, -1, 1,1);
	  manager -> AddProcess(new G4eIonisation,        -1, 2,2);
	  manager -> AddProcess(new G4eBremsstrahlung,    -1,-1,3);
	  manager -> AddProcess(new G4eplusAnnihilation,   0,-1,4); 
          manager -> AddProcess(new G4StepLimiter(),  -1,-1,3);
	}   
    }
}
