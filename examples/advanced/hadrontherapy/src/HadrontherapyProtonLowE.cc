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
// $Id: HadrontherapyProtonLowE.cc,v 1.1 2005-03-10 12:59:11 mpiergen Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "HadrontherapyProtonLowE.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4MultipleScattering.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4EnergyLossTables.hh"

HadrontherapyProtonLowE::HadrontherapyProtonLowE(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyProtonLowE::~HadrontherapyProtonLowE()
{ }

void HadrontherapyProtonLowE::ConstructProcess()
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
	       || particleName == "kaon-"

) 
	{

    G4hLowEnergyIonisation* ahadronLowEIon    = new G4hLowEnergyIonisation();
    ahadronLowEIon -> SetNuclearStoppingPowerModel("ICRU_R49") ; // ICRU49 models for nuclear SP
    ahadronLowEIon -> SetNuclearStoppingOn() ;
  
    // setting tables explicitly for electronic stopping power
    ahadronLowEIon -> SetElectronicStoppingPowerModel(G4GenericIon::GenericIonDefinition(), 
						    "ICRU_R49p") ;  // ICRU49 models for elettronic SP
    ahadronLowEIon -> SetElectronicStoppingPowerModel(G4Proton::ProtonDefinition(), 
						    "ICRU_R49p") ;
    // Switch off the Barkas and Bloch corrections
    ahadronLowEIon -> SetBarkasOff();

	    manager->AddProcess(new G4MultipleScattering,     -1,1,1); 
	    manager->AddProcess(new G4hLowEnergyIonisation,       -1,2,2);

	}   
    }
}
