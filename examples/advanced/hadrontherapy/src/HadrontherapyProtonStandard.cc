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
// $Id: HadrontherapyProtonStandard.cc,v 1.2 2005-04-28 20:39:33 mpiergen Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "HadrontherapyProtonStandard.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"

HadrontherapyProtonStandard::HadrontherapyProtonStandard(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyProtonStandard::~HadrontherapyProtonStandard()
{ }

void HadrontherapyProtonStandard::ConstructProcess()
{
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "proton"  
	       || particleName == "antiproton"  
	       || particleName == "pi+"  
	       || particleName == "pi-"  
	       || particleName == "kaon+"  
	       || particleName == "kaon-") 
	{

	    manager->AddProcess(new G4MultipleScattering,     -1,1,1); 
	    manager->AddProcess(new G4hIonisation,       -1,2,2);

	}   
    }
}
