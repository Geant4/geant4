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
// $Id: Tst50PhysicsList.cc,v 1.21 2004-07-20 15:11:12 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst50PhysicsList.hh"
#include "Tst50PhysicsListMessenger.hh"
#include "Tst50Particles.hh"
#include "Tst50AlphaICRU49.hh"
#include "Tst50AlphaStandard.hh"
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
#include "Tst50ProtonEEDL.hh"
#include "Tst50ProtonEEDLziegler.hh"
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

Tst50PhysicsList::Tst50PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false), 
                                      protonIsRegistered(false),
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
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50AlphaStandard(name) );
	  alphaIsRegistered = true;
	}
    }

//Proton
 if (name == "proton-eedl") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonEEDL(name) );
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

if (name == "proton-eedl-ziegler") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst50PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst50ProtonEEDLziegler(name) );
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










