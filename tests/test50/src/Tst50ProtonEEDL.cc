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
#include "Tst50ProtonEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4Proton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"

Tst50ProtonEEDL::Tst50ProtonEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50ProtonEEDL::~Tst50ProtonEEDL()
{ }

void Tst50ProtonEEDL::ConstructProcess()
{

theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "proton") 
	{
	  // G4VProcess*  multipleScattering= new G4MultipleScattering(); 
   manager->AddProcess(new G4hLowEnergyIonisation(),-1,2,2);
   //  manager->AddProcess(multipleScattering,-1,1,1);  	
    G4hLowEnergyLoss::SetEnlossFluc(false);
	}}}
