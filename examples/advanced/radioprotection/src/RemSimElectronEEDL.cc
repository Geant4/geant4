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
// $Id: RemSimElectronEEDL.cc,v 1.4 2005/05/19 13:46:29 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Author: Susanna Guatelli, guatelloi@ge.infn.it

#include "RemSimElectronEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"

RemSimElectronEEDL::RemSimElectronEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

RemSimElectronEEDL::~RemSimElectronEEDL()
{ }

void RemSimElectronEEDL::ConstructProcess()
{
  // Add EEDL processes for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if (particleName == "e-") 
	{
	  manager -> AddProcess(new G4MultipleScattering,     -1, 1,1);
	  manager -> AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	  manager -> AddProcess(new G4LowEnergyBremsstrahlung,-1,-1,3);
          manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}   
    }
}
