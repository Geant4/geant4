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
// $Id: RemSimElectronEEDL.cc,v 1.7 2009-11-12 05:12:18 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli, guatelloi@ge.infn.it

#include "RemSimElectronEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
// e-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
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
      G4ProcessManager* pmanager = particle -> GetProcessManager();
      G4String particleName = particle -> GetParticleName();
     
      if (particleName == "e-") 
	{
	  G4eMultipleScattering* msc = new G4eMultipleScattering();
	  msc->SetStepLimitType(fUseDistanceToBoundary);
	  pmanager->AddProcess(msc,-1, 1, 1);

	  // Ionisation
	  G4eIonisation* eIonisation = new G4eIonisation();
	  eIonisation->SetEmModel(new G4LivermoreIonisationModel());
	  eIonisation->SetStepFunction(0.2, 100*um); //improved precision in tracking  
	  pmanager->AddProcess(eIonisation,-1, 2, 2);
	
	  // Bremsstrahlung
	  G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
	  eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
	  pmanager->AddProcess(eBremsstrahlung, -1,-3, 3);

	  // Step Limiter
	  pmanager -> AddProcess(new G4StepLimiter(),  -1,-1,3);
	}   
    }
}
