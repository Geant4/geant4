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
// $Id: Tst50ElectronEEDLrange.cc,v 1.3 2003-07-28 15:05:52 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG          Designed for modular Physics List with
// CSDA and StoppingPower test conditions
//
// -------------------------------------------------------------------

#include "Tst50ElectronEEDLrange.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

Tst50ElectronEEDLrange::Tst50ElectronEEDLrange(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50ElectronEEDLrange::~Tst50ElectronEEDLrange()
{ }

void Tst50ElectronEEDLrange::ConstructProcess()
{
  // Add EEDLrangeprocesses for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "e-") 
	{
	  
	  manager->AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	  manager->AddProcess(new G4LowEnergyBremsstrahlung,-1,-1,3);
          G4VeLowEnergyLoss::SetEnlossFluc(false); 
	  G4cout<<" range conditions set: no msc, no energy fluct"<<G4endl;	
	}   
    }
}
