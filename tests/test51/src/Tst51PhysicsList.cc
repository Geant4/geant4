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
// $Id: Tst51PhysicsList.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
//
// -------------------------------------------------------------------

#include "Tst51PhysicsList.hh"
#include "Tst51PhysicsListMessenger.hh"
#include "Tst51Particles.hh"
#include "Tst51PhotonStandard.hh"
#include "Tst51PhotonEPDL.hh"
#include "Tst51PhotonPenelope.hh"
#include "Tst51PhotonPolarised.hh"
#include "Tst51ElectronStandard.hh"
#include "Tst51ElectronEEDL.hh"
#include "Tst51ElectronEEDL2BS.hh"
#include "Tst51ElectronEEDL2BN.hh"
#include "Tst51ElectronPenelope.hh"
#include "Tst51PositronStandard.hh"
#include "Tst51PositronPenelope.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"


Tst51PhysicsList::Tst51PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false)                                      
{
  defaultCutValue = 0.5 * micrometer;
  SetVerboseLevel(1);

  messenger = new Tst51PhysicsListMessenger(this);
 
  RegisterPhysics( new Tst51Particles("particles") );
}


Tst51PhysicsList::~Tst51PhysicsList()
{
  delete messenger;
}


void Tst51PhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51PhotonStandard(name) );
	  photonIsRegistered = true;
	}
    }
  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51PhotonEPDL(name) );
	  photonIsRegistered = true;
	}
   } 
  // Register processes a' la Penelope for photons
  if (name == "photon-penelope")
    {
     if (photonIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51PhotonPenelope(name) );
	   photonIsRegistered = true;
	}
    }

  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51ElectronStandard(name) );	  
	  electronIsRegistered = true;
	}
    }

  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl-tsai") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51ElectronEEDL(name) );
	  electronIsRegistered = true;
	}
   } 

if (name == "electron-eedl-2bn") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51ElectronEEDL2BN(name) );
	  electronIsRegistered = true;
	}
   } 


if (name == "electron-eedl-2bs") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51ElectronEEDL2BS(name) );
	  electronIsRegistered = true;
	}
   } 

  // Register processes a' la Penelope for electrons
  if (name == "electron-penelope")
    {
     if (electronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51ElectronPenelope(name) );
	  electronIsRegistered = true;
	}
    }
  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51PositronStandard(name) );
	  positronIsRegistered = true;
	}
    }

 if (name == "positron-penelope")
    {
     if (positronIsRegistered) 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst51PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst51PositronPenelope(name) );
	  positronIsRegistered = true;
	}
    }

  if (electronIsRegistered && positronIsRegistered && photonIsRegistered)
    {
      G4cout << "PhysicsList for electron, positron and photon registered" << G4endl;
    }
}

void Tst51PhysicsList::SetParticleCut(G4double value)
{
  defaultCutValue = value; 
  G4cout << "The threshold of production of secondaries is now (cm):" 
         << value/cm << G4endl; 
}

void Tst51PhysicsList::SetCuts()
{
 G4VUserPhysicsList::SetCutsWithDefault();
 if (verboseLevel>0) DumpCutValuesTable();
}










