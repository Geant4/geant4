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
// $Id: Tst14ElectronEEDL.cc,v 1.5 2010-09-01 17:29:29 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria.Grazia.Pia@cern.ch
//
// History:
// -----------
// 22 Feb 2003 MGP          Designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst14ElectronEEDL.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

Tst14ElectronEEDL::Tst14ElectronEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst14ElectronEEDL::~Tst14ElectronEEDL()
{ }

void Tst14ElectronEEDL::ConstructProcess()
{
  // Add EEDL processes for electrons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      G4double LivermoreHighEnergyLimit = GeV;

      if (particleName == "e-") 
	{

	  G4eMultipleScattering* msc = new G4eMultipleScattering();
	  //msc->AddEmModel(0, new G4UrbanMscModel93());
	  msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
	  msc->SetStepLimitType(fUseDistanceToBoundary);
	  manager->AddProcess(msc,                   -1, 1, 1);
      
	  // Ionisation
	  G4eIonisation* eIoni = new G4eIonisation();
	  G4LivermoreIonisationModel* theIoniLivermore = new
	    G4LivermoreIonisationModel();
	  theIoniLivermore->SetHighEnergyLimit(1*MeV); 
	  eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
	  eIoni->SetStepFunction(0.2, 100*um); //     
	  manager->AddProcess(eIoni,                 -1, 2, 2);
      
	  // Bremsstrahlung
	  G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
	  G4LivermoreBremsstrahlungModel* theBremLivermore = new
	    G4LivermoreBremsstrahlungModel();
	  theBremLivermore->SetHighEnergyLimit(LivermoreHighEnergyLimit);
	  eBrem->AddEmModel(0, theBremLivermore);
	  manager->AddProcess(eBrem, -1,-3, 3);
	
	}   
    }
}
