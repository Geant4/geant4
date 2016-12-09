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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPhysicsList.cc                           *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models 
// 
//
// **********************************************************************

#include "XrayTelPhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
// Hadrons
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

XrayTelPhysicsList::XrayTelPhysicsList():  G4VUserPhysicsList()
{
  // Default cut values
  defaultCutValue = 2.0*mm;
  cutForGamma     = 1.0*m;
  cutForElectron  = 1.0*m;

  SetVerboseLevel(1);
}

XrayTelPhysicsList::~XrayTelPhysicsList()
{}

void XrayTelPhysicsList::ConstructParticle()
{
  // Here are constructed all particles
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
  ConstructAllShortLiveds();
}

// In this method, static member functions should be called for ALL particles to be used.

void XrayTelPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
void XrayTelPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  G4MuonPlus::Definition();
  G4MuonMinus::Definition();
}
void XrayTelPhysicsList::ConstructMesons()
{
  // mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
}

void XrayTelPhysicsList::ConstructBaryons()
{
  //  barions
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();
}

void XrayTelPhysicsList::ConstructIons()
{
  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

void XrayTelPhysicsList::ConstructAllShortLiveds()
{
  //Short-lived
  G4ShortLivedConstructor slConstructor;
  slConstructor.ConstructParticle();
}


void XrayTelPhysicsList::ConstructProcess()
{
  // Transportation, electromagnetic and general processes 
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

// Here are respective header files for chosen processes

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

void XrayTelPhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();

  while( (*particleIterator)() )
    {
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma")
	{
	  //gamma
	  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
	  pmanager->AddDiscreteProcess(new G4ComptonScattering());
	  pmanager->AddDiscreteProcess(new G4GammaConversion());
	}
      else if (particleName == "e-")
	{
	  //electron
	  pmanager->AddProcess(new G4eMultipleScattering(),-1, 1,1);
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3,3);

	} else if (particleName == "e+") {

	  //positron
	  pmanager->AddProcess(new G4eMultipleScattering(),-1, 1,1);
	  pmanager->AddProcess(new G4eIonisation(),       -1, 2,2);
	  pmanager->AddProcess(new G4eBremsstrahlung(),   -1, 3,3);
	  pmanager->AddProcess(new G4eplusAnnihilation(),  0,-1,4);

	} else if( particleName == "mu-" || 
		   particleName == "mu+" ) {

	  //muon  
	  pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);
	  pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 3, 3);
	  pmanager->AddProcess(new G4MuPairProduction,    -1, 4, 4);
	  
	} else if( particleName == "pi-" ||
		   particleName == "pi+" ) {

	  //pions
	  pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
	  pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
	  pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);       
     
	} else if( particleName == "proton" ) {

	  //proton  
	  pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
	  pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
	  pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);       
     
	} else if( particleName == "alpha" || 
		   particleName == "He3" || 
		   particleName == "GenericIon" ) {

	  //Ions 
	  pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
	  
	} else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino")) {

	  //all others charged particles except geantino
	  pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
	  pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);

	}
    }
}

#include "G4Decay.hh"

void XrayTelPhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void XrayTelPhysicsList::SetCuts()
{
  // defaultCutValue you have typed in is used

  if (verboseLevel >1){
    G4cout << "XrayTelPhysicsList::SetCuts:" << G4endl;
  }

  // set cut values for gamma at first and for e- second
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForElectron, "e+");

  // SetCutValueForOthers(defaultCutValue); 
 
  if (verboseLevel >1) { 
    DumpCutValuesTable(); 
  } 
}

void XrayTelPhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

void XrayTelPhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

G4double XrayTelPhysicsList::GetCutForGamma() const
{
  return cutForGamma;
}

G4double XrayTelPhysicsList::GetCutForElectron() const
{
  return cutForElectron;
}





