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
// $Id: HadrontherapyPhysicsList.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "globals.hh"
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
//#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh" 
//#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4VeLowEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include <iomanip.h>               


HadrontherapyPhysicsList::HadrontherapyPhysicsList(): G4VModularPhysicsList(),
						      
electronIsRegistered(false), positronIsRegistered(false), 
photonIsRegistered(false), ionIsRegistered(false),
protonPrecompoundIsRegistered(false), muonIsRegistered(false),
decayIsRegistered(false)
{
  defaultCutValue = 10. * mm;
  currentDefaultCut = defaultCutValue;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton = defaultCutValue;
  messenger = new HadrontherapyPhysicsListMessenger(this);
  SetVerboseLevel(1);
}

HadrontherapyPhysicsList::~HadrontherapyPhysicsList()
{ 
  delete messenger;
}

void HadrontherapyPhysicsList::AddPhysicsList(const G4String& name)
{
  G4cout << "Adding PhysicsList chunk " << name << G4endl;

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
		 << " cannot be registered ---- electron List already existing"                  << G4endl;
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
		 << " cannot be registered ---- electron List already existing"                  << G4endl;
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
		 << " cannot be registered ---- electron List already existing"                  << G4endl;
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
		 << " cannot be registered ---- positron List already existing"                  << G4endl;
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

  // Register Low Energy ICRU processes for charged particles
  //excluding e+, e-, muons
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

 // Register Standard processes for charged particles
  //excluding e+, e-, muons

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

  if (electronIsRegistered && positronIsRegistered && photonIsRegistered &&
      ionIsRegistered) 
    {
     G4cout<< 
     "Electromagnetic physics is registered for electron, positron, photons, protons" 
     << G4endl;
    }
    if (protonPrecompoundIsRegistered && muonIsRegistered && decayIsRegistered)
      {
	G4cout << " Hadronic physics is registered"<< G4endl;
      }     
}

void HadrontherapyPhysicsList::SetCuts()
{ 
 
 // reactualise cutValues
  if (currentDefaultCut != defaultCutValue)
    {
     if(cutForGamma    == currentDefaultCut) cutForGamma    = defaultCutValue;
     if(cutForElectron == currentDefaultCut) cutForElectron = defaultCutValue;
     if(cutForProton   == currentDefaultCut) cutForProton   = defaultCutValue;
     currentDefaultCut = defaultCutValue;
    }
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 

  //SetCutValue(defaultCutValue, "proton");
  //SetCutValue(defaultCutValue, "anti_proton");

  
  // Cut per region
  // in region dosemeter we need a very accurate precision

  G4String regName = "PhantomLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regName);
  G4ProductionCuts* cuts = new G4ProductionCuts ;
  cuts -> SetProductionCut(0.001*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(0.001*mm,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(0.001*mm,G4ProductionCuts::GetIndex("e+"));
  region -> SetProductionCuts(cuts);
      
  //if (verboseLevel >0){
  //    G4cout << "HadrontherapyPhysicsList::SetCuts:";
  //    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") 
  //    << G4endl;
  //  }
    
    if (verboseLevel>0) DumpCutValuesTable();
}

void HadrontherapyPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

void HadrontherapyPhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

void HadrontherapyPhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

// ---------------------------------------------------------------------------
//void HadrontherapyPhysicsList::GetRange(G4double val)
//{
//G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
//G4Material* currMat = pDet->GetPhantomMaterial();
//
//G4ParticleDefinition* part;
//G4double cut;
//part = theParticleTable->FindParticle("e-");
//cut = G4EnergyLossTables::GetRange(part,val,currMat);
//G4cout << "material : " << currMat->GetName() << G4endl;
//G4cout << "particle : " << part->GetParticleName() << G4endl;
//G4cout << "energy   : " << G4BestUnit(val,"Energy") << G4endl;
//G4cout << "range    : " << G4BestUnit(cut,"Length") << G4endl;
//}


