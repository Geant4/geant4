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
// $Id: HadrontherapyPhysicsList.cc,v 1.0
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyParticles.hh"
#include "HadrontherapyPhotonStandard.hh"
#include "HadrontherapyPhotonEPDL.hh"
#include "HadrontherapyPhotonPenelope.hh"
#include "HadrontherapyElectronStandard.hh"
#include "HadrontherapyElectronEEDL.hh"
#include "HadrontherapyElectronPenelope.hh"
#include "HadrontherapyPositronStandard.hh"
#include "HadrontherapyPositronPenelope.hh"
#include "HadrontherapyIonLowE.hh"
#include "HadrontherapyIonStandard.hh"
#include "HadrontherapyProtonPrecompound.hh"
#include "HadrontherapyMuonStandard.hh"
#include "HadrontherapyDecay.hh"

HadrontherapyPhysicsList::HadrontherapyPhysicsList(): G4VModularPhysicsList(),
						      electronIsRegistered(false), 
						      positronIsRegistered(false), 
						      photonIsRegistered(false), 
						      ionIsRegistered(false),
						      protonPrecompoundIsRegistered(false), 
						      muonIsRegistered(false),
						      decayIsRegistered(false)
{
  // The threshold of production of secondaries is fixed to 10. mm
  // for all the particles, in all the experimental set-up
  defaultCutValue = 10. * mm;
  messenger = new HadrontherapyPhysicsListMessenger(this);
  SetVerboseLevel(1);
}

HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{ 
  delete messenger;
}

void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{
  G4cout << "Adding PhysicsList component " << name << G4endl;
  
  //
  // Electromagnetic physics. 
  //
  // The user can choose three alternative approaches:
  // Standard, Low Energy based on the Livermore libraries and Low Energy Penelope
  //

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyPhotonStandard(name) );
	  photonIsRegistered = true;
	}
    }

  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing"
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyPhotonEPDL(name) );
	  photonIsRegistered = true;
	}
    } 

  // Register processes a' la Penelope for photons
  if (name == "photon-penelope")
    {
      if (photonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyPhotonPenelope(name) );
	  photonIsRegistered = true;
	}
    }

  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" 
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyElectronStandard(name) );	  
	  electronIsRegistered = true;
	}
    }

  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyElectronEEDL(name) );
	  electronIsRegistered = true;
	}
    } 

  // Register processes a' la Penelope for electrons
  if (name == "electron-penelope")
    {
      if (electronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyElectronPenelope(name) );
	  electronIsRegistered = true;
	}
    }

  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing"                  
		 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyPositronStandard(name) );
	  positronIsRegistered = true;
	}
    }
  // Register penelope processes for positrons
  if (name == "positron-penelope") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing"                  << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyPositronPenelope(name) );
	  positronIsRegistered = true;
	}
    }

  // Register Low Energy ICRU processes for protons and ions

  if (name == "ion-LowE") 
    {
      if (ionIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyIonLowE(name) );
	  ionIsRegistered = true;
	}
    }

  // Register Standard processes for protons and ions

  if (name == "ion-standard") 
    {
      if (ionIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- ion List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyIonStandard(name) );
	  ionIsRegistered = true;
	}
    }

  // Register the Standard processes for muons
  if (name == "muon-standard") 
    {
      if (muonIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyMuonStandard(name) );
	  muonIsRegistered = true;
	}
    }

  // Register the decay process
  if (name == "decay") 
    {
      if (decayIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyDecay(name) );
	  decayIsRegistered = true;
	}
    }
  //
  //
  // Register the hadronic physics for protons, neutrons, ions
  //
  // 
  if (name == "proton-precompound") 
    {
      if (protonPrecompoundIsRegistered) 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "HadrontherapyPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new HadrontherapyProtonPrecompound(name) );
	  protonPrecompoundIsRegistered = true;
	}
    }

  if (electronIsRegistered && positronIsRegistered && photonIsRegistered &&
      ionIsRegistered) 
    {
      G4cout << 
	"Electromagnetic physics is registered for electron, positron, photons, protons" 
	     << G4endl;
    }
    if (protonPrecompoundIsRegistered && muonIsRegistered && decayIsRegistered)
      {
	G4cout << " Hadronic physics is registered" << G4endl;
      }     
}

void HadrontherapyPhysicsList::SetCuts()
{  
  // Set the threshold of production equal to the defaultCutValue
  // in the experimental set-up
  G4VUserPhysicsList::SetCutsWithDefault();
    
  // Definition of a smaller threshold of production in the phantom region
  // where high accuracy is required in the energy deposit calculation

  G4String regionName = "PhantomLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);
  G4ProductionCuts* cuts = new G4ProductionCuts ;
  G4double regionCut = 0.001*mm;
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e+"));
  region -> SetProductionCuts(cuts);

  if (verboseLevel>0) DumpCutValuesTable();
}


