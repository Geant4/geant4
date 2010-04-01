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
// Author: Stephane Chauvie (chauvie@to.infn.it)
//
// History:
// -----------
// 22 Feb 2006         SC      Created
//
// -------------------------------------------------------------------
                                                                                                            
#include "Tst50AntiProtonICRU49.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
//#include "G4MultipleScattering.hh"
#include "G4AntiProton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4StepLimiter.hh"

Tst50AntiProtonICRU49::Tst50AntiProtonICRU49(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50AntiProtonICRU49::~Tst50AntiProtonICRU49()
{ }

void Tst50AntiProtonICRU49::ConstructProcess()
{

theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "anti_proton" )
	{
	  G4hLowEnergyIonisation* ionisation = new G4hLowEnergyIonisation();
          // G4VProcess*  multipleScattering= new G4MultipleScattering(); 
	  ionisation -> SetEnlossFluc(false); 

          ionisation -> SetNuclearStoppingOn() ;
          ionisation -> SetBarkasOn();

          manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	  manager -> AddProcess(ionisation,-1,2,2);
          //  manager->AddProcess(multipleScattering,-1,1,1);  	
	}	
    }
}
