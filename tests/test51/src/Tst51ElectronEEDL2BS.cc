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
// $Id: Tst51ElectronEEDL2BS.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst51ElectronEEDL2BS.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"
#include "G4VBremAngularDistribution.hh"
#include "G4Generator2BS.hh"

Tst51ElectronEEDL2BS::Tst51ElectronEEDL2BS(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst51ElectronEEDL2BS::~Tst51ElectronEEDL2BS()
{ }

void Tst51ElectronEEDL2BS::ConstructProcess()
{
  // Add EEDL  processes for electrons with 2BS angular distribution
  
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
          
	  */
           manager->AddProcess(new G4MultipleScattering,     -1, 1,1);
	   manager->AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	   G4LowEnergyBremsstrahlung* brem = new G4LowEnergyBremsstrahlung();
           brem -> SetAngularGenerator("2bs");  
           manager -> AddProcess(brem,-1,-1,3);
	}   
    }
}
