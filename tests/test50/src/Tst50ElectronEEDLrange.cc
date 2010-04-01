//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: Tst50ElectronEEDLrange.cc,v 1.7 2010-04-01 09:48:30 sincerti Exp $
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
#include "G4eMultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"

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
	  
	  G4LowEnergyIonisation* loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
          G4LowEnergyBremsstrahlung* loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
          loweIon -> SetEnlossFluc(false); 
          loweBrem -> SetEnlossFluc(false);
          manager -> AddProcess(loweIon,     -1, 2,2);
          manager -> AddProcess(loweBrem,    -1,-1,3);  
          manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
          G4cout<<" range conditions set: no msc, no energy fluct"<<G4endl;	
	}   
    }
}
