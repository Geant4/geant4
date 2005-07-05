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
// $Id: Tst51ElectronEEDL.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst51ElectronEEDL.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"
#include "G4VBremAngularDistribution.hh"

Tst51ElectronEEDL::Tst51ElectronEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst51ElectronEEDL::~Tst51ElectronEEDL()
{ }

void Tst51ElectronEEDL::ConstructProcess()
{
  // Add EEDL processes for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "e-") 
	{
	  /* physics list for premiliminary study
	  // manager->AddProcess(new G4MultipleScattering,     -1, 1,1);
	  //manager->AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	  G4LowEnergyBremsstrahlung* brem = new G4LowEnergyBremsstrahlung();
	  brem -> SetEnlossFluc(false); 
          //brem -> SetAngularGenerator("2bn");  
              
	   brem -> SetAngularGenerator("2bs");  
          manager->AddProcess(brem,-1,-1,3);
          
	  //manager->AddProcess(new G4StepLimiter(),          -1,-1, 3);
	  */
           manager->AddProcess(new G4MultipleScattering,     -1, 1,1);
	   manager->AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	   G4LowEnergyBremsstrahlung* brem = new G4LowEnergyBremsstrahlung();
	 	   
           brem -> SetAngularGenerator("tsai");  
           manager -> AddProcess(brem,-1,-1,3);
	}   
    }
}
