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
// $Id: Tst50PhysicsList.cc,v 1.26 2010-06-25 09:46:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
// 22 Feb 2005 SC           Added Antiproton Processes
//
// -------------------------------------------------------------------

#include "Tst50PhysicsList.hh"
#include "Tst50PhysicsListMessenger.hh"
#include "Tst50Particles.hh"
#include "Tst50AlphaICRU49.hh"
#include "Tst50AlphaStandard.hh"
#include "Tst50AlphaZiegler.hh"
#include "Tst50PhotonStandard.hh"
#include "Tst50PhotonEPDL.hh"
#include "Tst50PhotonPenelope.hh"
#include "Tst50PhotonPolarised.hh"
#include "Tst50ElectronStandard.hh"
#include "Tst50ElectronStandardback.hh"
#include "Tst50ElectronEEDL.hh"
#include "Tst50ElectronEEDLrange.hh"
#include "Tst50ElectronEEDLback.hh"
#include "Tst50ElectronPenelope.hh"
#include "Tst50PositronStandard.hh"
#include "Tst50PositronStandardBack.hh"
#include "Tst50PositronPenelope.hh"
#include "Tst50ProtonStandard.hh"
#include "Tst50ProtonICRU49.hh"
#include "Tst50ProtonZiegler85.hh"
#include "Tst50AntiProtonICRU49.hh"
#include "Tst50AntiProtonZiegler85.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "Tst50ProtonZiegler2000.hh"
#include "Tst50AntiProtonZiegler2000.hh"

Tst50PhysicsList::Tst50PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false), 
                                      protonIsRegistered(false),
                                      anti_protonIsRegistered(false),
                                      alphaIsRegistered(false)
{
  defaultCutValue = 0.1 * mm;
  SetVerboseLevel(1);

  messenger = new Tst50PhysicsListMessenger(this);
 
  RegisterPhysics( new Tst50Particles("particles") );
}


Tst50PhysicsList::~Tst50PhysicsList()
{
  delete messenger;
}


void Tst50PhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PhotonStandard(name) );
	  photonIsRegistered = true;
	}
    }
  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PhotonEPDL(name) );
	  photonIsRegistered = true;
	}
   } 
  // Register processes a' la Penelope for photons
  if (name == "photon-penelope")
    {
     if (photonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PhotonPenelope(name) );
	   photonIsRegistered = true;
	}
    }

  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronStandard(name) );	  
	  electronIsRegistered = true;
	}
    }
 if (name == "electron-standard-back") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronStandardback(name) );	  
	  electronIsRegistered = true;
	}
    }
  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronEEDL(name) );
	  electronIsRegistered = true;
	}
   } 
 if (name == "electron-eedl-back") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronEEDLback(name) );
	  electronIsRegistered = true;
	}
   } 

 if (name == "electron-eedl-range") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronEEDLrange(name) );
	  electronIsRegistered = true;
	}
   } 
  // Register processes a' la Penelope for electrons
  if (name == "electron-penelope")
    {
     if (electronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ElectronPenelope(name) );
	  electronIsRegistered = true;
	}
    }
  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PositronStandard(name) );
	  positronIsRegistered = true;
	}
    }


 if (name == "positron-standard-back") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PositronStandardBack(name) );
	  protonIsRegistered = true;
	}
    }

 if (name == "positron-penelope")
    {
     if (positronIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50PositronPenelope(name) );
	  positronIsRegistered = true;
	}
    }

 //Alpha particle  
if (name == "alpha-ICRU49") 
    {
      if (alphaIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- alpha e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AlphaICRU49(name) );
	  alphaIsRegistered = true;
	}
    }

if (name == "alpha-standard") 
    {
      if (alphaIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- alpha e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name
           << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AlphaStandard(name) );
	  alphaIsRegistered = true;
	}
    }

if (name == "alpha-ziegler77") 
    {
      if (alphaIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- alpha e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AlphaZiegler(name) );
	  alphaIsRegistered = true;
	}
    }

// Proton

 if (name == "proton-ICRU49") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonICRU49(name) );
	  protonIsRegistered = true;
	}
    }
 if (name == "proton-ziegler2000") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonZiegler2000(name) );
	  protonIsRegistered = true;
	}
    }

if (name == "proton-ziegler85") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonZiegler85(name) );
	  protonIsRegistered = true;
	}
    }

if (name == "proton-standard") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonStandard(name) );
	  protonIsRegistered = true;
	}
    }
// Anti Proton

 if (name == "anti_proton-ICRU49") 
    {
      if (anti_protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- antiproton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AntiProtonICRU49(name) );
	  anti_protonIsRegistered = true;
	}
    }
 if (name == "anti_proton-ziegler2000") 
    {
      if (anti_protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- anti proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AntiProtonZiegler2000(name) );
	  anti_protonIsRegistered = true;
	}
    }

if (name == "anti_proton-ziegler85") 
    {
      if (anti_protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- anti proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AntiProtonZiegler85(name) );
	  anti_protonIsRegistered = true;
	}
    }


  if (electronIsRegistered && positronIsRegistered && photonIsRegistered)
    {
      G4cout << "PhysicsList for electron, positron and photon registered" << G4endl;
    }
}

void Tst50PhysicsList::SetGammaCut(G4double value)
{
  ResetCuts();
  cutForGamma = value;
}


void Tst50PhysicsList::SetElectronCut(G4double value)
{
  ResetCuts();
  cutForElectron = value;
}

void Tst50PhysicsList::SetParticleCut(G4double value)
{
  defaultCutValue = value; 
}

void Tst50PhysicsList::SetCuts()
{
 G4VUserPhysicsList::SetCutsWithDefault();
 if (verboseLevel>0) DumpCutValuesTable();
}










