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
// $Id: Tst50ProtonEEDLziegler.cc,v 1.4 2004-07-20 15:11:12 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG          Designed for modular Physics List with
// Ziegler model for the Stopping Power and CSDArange and Stopping Power 
//conditions set
//
// -------------------------------------------------------------------

#include "Tst50ProtonEEDLziegler.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4Proton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4hZiegler1985p.hh"

Tst50ProtonEEDLziegler::Tst50ProtonEEDLziegler(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50ProtonEEDLziegler::~Tst50ProtonEEDLziegler()
{ }

void Tst50ProtonEEDLziegler::ConstructProcess()
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
	  G4hLowEnergyIonisation* ion= new G4hLowEnergyIonisation();
	  
	  ion -> SetElectronicStoppingPowerModel(particle, "Ziegler1985p");
          ion -> SetNuclearStoppingPowerModel("Ziegler1985");
          ion -> SetNuclearStoppingOn() ;
	  //  manager->AddProcess(multipleScattering,-1,1,1);  	
	  ion -> SetEnlossFluc(false);
          
          manager->AddProcess(ion,-1,2,2);
	}
    }
}
