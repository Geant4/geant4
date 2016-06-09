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
//
// Code developed by:
// S. Guatelli, G. Guerrieri and M. G. Pia
//
// 
//
//    **********************************
//    *                                *
//    *     G4HumanPhantomPhysicsList.cc       *
//    *                                *
//    **********************************
//
// $Id$
//
#include "G4HumanPhantomPhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"              
// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

// e-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

G4HumanPhantomPhysicsList::G4HumanPhantomPhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

G4HumanPhantomPhysicsList::~G4HumanPhantomPhysicsList()
{
}

void G4HumanPhantomPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();

}

void G4HumanPhantomPhysicsList::ConstructBosons()
{ 
  // photons
  G4Gamma::GammaDefinition();
}

void G4HumanPhantomPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void G4HumanPhantomPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

void G4HumanPhantomPhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    // Processes
    
    if (particleName == "gamma") {
      
      // Photon     
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      theRayleigh->SetModel(new G4LivermoreRayleighModel());  //not strictly necessary
      pmanager->AddDiscreteProcess(theRayleigh);

      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetModel(new G4LivermorePhotoElectricModel());
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
	
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetModel(new G4LivermoreComptonModel());
      pmanager->AddDiscreteProcess(theComptonScattering);
	
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetModel(new G4LivermoreGammaConversionModel());
      pmanager->AddDiscreteProcess(theGammaConversion);
      
    } else if (particleName == "e-") {
      // Electron
      
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
      
    } else if (particleName == "e+") {
      // Positron      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,-1, 1, 1);
	
      // Ionisation
      G4eIonisation* eIonisation = new G4eIonisation();
      eIonisation->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIonisation,                 -1, 2, 2);

      //Bremsstrahlung (use default, no low-energy available)
      pmanager->AddProcess(new G4eBremsstrahlung(), -1,-1, 3);

      //Annihilation
      pmanager->AddProcess(new G4eplusAnnihilation(),0,-1, 4);     
      
    }
  }  
}

void G4HumanPhantomPhysicsList::SetCuts()
{
  // The production threshold is fixed to 0.1 mm for all the particles
  // Secondary particles with a range bigger than 0.1 mm 
  // are generated; otherwise their energy is considered deposited locally

  defaultCutValue = 0.1 * mm;

  const G4double cutForGamma = defaultCutValue;
  const G4double cutForElectron = defaultCutValue;
  const G4double cutForPositron = defaultCutValue;

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  // Set the secondary production cut lower than 990. eV
  // Very important for high precision of lowenergy processes at low energies
 
  G4double lowLimit = 250. * eV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
  
  if (verboseLevel>0) DumpCutValuesTable();
}
