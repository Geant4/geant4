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
// $Id: RemSimPhysicsList.cc,v 1.2 2004-02-03 09:16:47 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
//
// -------------------------------------------------------------------

#include "RemSimPhysicsList.hh"
#include "RemSimPhysicsListMessenger.hh"
#include "RemSimParticles.hh"
#include "RemSimPhotonStandard.hh"
#include "RemSimPhotonEPDL.hh"
#include "RemSimElectronStandard.hh"
#include "RemSimElectronEEDL.hh"
#include "RemSimPositronStandard.hh"
#include "RemSimProtonStandard.hh"
#include "RemSimProtonEEDL.hh"
#include "RemSimProtonEEDLziegler.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

RemSimPhysicsList::RemSimPhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false), 
                                      protonIsRegistered(false)
{
  defaultCutValue = 0.1 * mm;
  SetVerboseLevel(1);

  messenger = new RemSimPhysicsListMessenger(this);
 
  RegisterPhysics( new RemSimParticles("particles") );
}


RemSimPhysicsList::~RemSimPhysicsList()
{
  delete messenger;
}


void RemSimPhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimPhotonStandard(name) );
	  photonIsRegistered = true;
	}
    }
  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimPhotonEPDL(name) );
	  photonIsRegistered = true;
	}
    } 

  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimElectronStandard(name) );	  
	  electronIsRegistered = true;
	}
    }
  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimElectronEEDL(name) );
	  electronIsRegistered = true;
	}
   } 

  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimPositronStandard(name) );
	  positronIsRegistered = true;
	}
    }

 if (name == "proton-eedl") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimProtonEEDL(name) );
	  protonIsRegistered = true;
	}
    }

if (name == "proton-eedl-ziegler") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimProtonEEDLziegler(name) );
	  protonIsRegistered = true;
	}
    }

if (name == "proton-standard") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new RemSimProtonStandard(name) );
	  protonIsRegistered = true;
	}
    }

  if (electronIsRegistered && positronIsRegistered && photonIsRegistered)
    {
      G4cout << "PhysicsList for electron, positron and photon registered" << G4endl;
    }
}

void RemSimPhysicsList::SetGammaCut(G4double value)
{
  ResetCuts();
  cutForGamma = value;
}


void RemSimPhysicsList::SetElectronCut(G4double value)
{
  ResetCuts();
  cutForElectron = value;
}

void RemSimPhysicsList::SetParticleCut(G4double value)
{
  defaultCutValue = value; 
}

void RemSimPhysicsList::SetCuts()
{
 G4VUserPhysicsList::SetCutsWithDefault();
 if (verboseLevel>0) DumpCutValuesTable();
}










