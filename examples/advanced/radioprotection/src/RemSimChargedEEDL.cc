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
// $Id: RemSimChargedEEDL.cc,v 1.1 2004-03-12 10:56:53 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "RemSimChargedEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4EnergyLossTables.hh"
#include "G4MultipleScattering.hh"

RemSimChargedEEDL::RemSimChargedEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

RemSimChargedEEDL::~RemSimChargedEEDL()
{ }

void RemSimChargedEEDL::ConstructProcess()
{
  // Add standard processes for muons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      G4double charge = particle->GetPDGCharge();

      if(particleName != "proton" && particleName != "alpha" &&
         particleName != "e+" &&  particleName != "e-")
	{ 
      if(charge!= 0.0) 
      {
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
        G4hLowEnergyIonisation* ahadronLowEIon = new G4hLowEnergyIonisation();
	pmanager->AddProcess(aMultipleScattering,-1,1,1);
	pmanager->AddProcess(ahadronLowEIon,     -1,2,2);   
      } 
	}
    }
}
