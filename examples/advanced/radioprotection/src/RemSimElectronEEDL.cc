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
// $Id: RemSimElectronEEDL.cc,v 1.6 2006/06/29 16:23:43 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Author: Susanna Guatelli, guatelloi@ge.infn.it

#include "RemSimElectronEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"

RemSimElectronEEDL::RemSimElectronEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

RemSimElectronEEDL::~RemSimElectronEEDL()
{ }

void RemSimElectronEEDL::ConstructProcess()
{
  // Add EEDL processes for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator -> value();
      G4ProcessManager* manager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if (particleName == "e-") 
	{
	  manager -> AddProcess(new G4MultipleScattering,     -1, 1,1);
	  manager -> AddProcess(new G4LowEnergyIonisation,    -1, 2,2);
	  manager -> AddProcess(new G4LowEnergyBremsstrahlung,-1,-1,3);
          manager -> AddProcess(new G4StepLimiter(),-1,-1,3);
	}   
    }
}
