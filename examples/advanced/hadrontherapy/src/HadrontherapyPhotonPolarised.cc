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
// $Id: HadrontherapyPhotonPolarised.cc,v 1.1 2005-03-10 12:59:11 mpiergen Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// -------------------------------------------------------------------

#include "HadrontherapyPhotonPolarised.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

HadrontherapyPhotonPolarised::HadrontherapyPhotonPolarised(const G4String& name): G4VPhysicsConstructor(name)
{ }

HadrontherapyPhotonPolarised::~HadrontherapyPhotonPolarised()
{ }

void HadrontherapyPhotonPolarised::ConstructProcess()
{
  // Add polarised processes for photons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "gamma") 
	{
	  manager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
	  manager->AddDiscreteProcess(new G4LowEnergyPolarizedCompton);
	  manager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
	  manager->AddDiscreteProcess(new G4LowEnergyRayleigh);
	}   
    }
}
