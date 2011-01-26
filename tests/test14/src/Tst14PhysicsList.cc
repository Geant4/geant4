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
// $Id: Tst14PhysicsList.cc,v 1.25 2010-08-29 19:50:17 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
// 15 Dec 2008 L. Pandola   Added constructor for e+ Penelope
//
// -------------------------------------------------------------------

#include "Tst14PhysicsList.hh"
#include "Tst14PhysicsListMessenger.hh"
#include "Tst14Particles.hh"
#include "Tst14PhotonStandard.hh"
#include "Tst14PhotonEPDL.hh"
#include "Tst14PhotonPenelope.hh"
#include "Tst14PhotonPolarised.hh"
#include "Tst14ElectronStandard.hh"
#include "Tst14ElectronEEDL.hh"
#include "Tst14ElectronPenelope.hh"
#include "Tst14PositronStandard.hh"
#include "Tst14PositronPenelope.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"



Tst14PhysicsList::Tst14PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false)
{
  defaultCutValue = 0.01 * um;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;

  SetVerboseLevel(1);

  // UI messenger
  messenger = new Tst14PhysicsListMessenger(this);
 
  // Particles
  RegisterPhysics( new Tst14Particles("particles") );

}


Tst14PhysicsList::~Tst14PhysicsList()
{
  delete messenger;
}


void Tst14PhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PhotonStandard(name) );
	  photonIsRegistered = true;
	}
    }
  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PhotonEPDL(name) );
	  photonIsRegistered = true;
	}
   } 
  // Register processes a' la Penelope for photons
  if (name == "photon-penelope")
    {
     if (photonIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PhotonPenelope(name) );
	  photonIsRegistered = true;
	}
    }
  // Register polarised processes for photons
  if (name == "photon-polarised")
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PhotonPolarised(name) );
	  photonIsRegistered = true;
	}
    }
  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14ElectronStandard(name) );	  
	  electronIsRegistered = true;
	}
    }
  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14ElectronEEDL(name) );
	  electronIsRegistered = true;
	}
   } 
  // Register processes a' la Penelope for electrons
  if (name == "electron-penelope")
    {
     if (electronIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14ElectronPenelope(name) );
	  electronIsRegistered = true;
	}
    }
  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PositronStandard(name) );
	  positronIsRegistered = true;
	}
    }

  if (name == "positron-penelope") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst14PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst14PositronPenelope(name) );
	  positronIsRegistered = true;
	}
    }


  if (electronIsRegistered && positronIsRegistered && photonIsRegistered)
    {
      G4cout << "PhysicsList for electron, positron and photon registered" << G4endl;
    }
}


void Tst14PhysicsList::SetGELowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Gamma and electron cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(cut,1e5);
}


void Tst14PhysicsList::SetGammaLowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Gamma cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(cut,1e5);
}

void Tst14PhysicsList::SetElectronLowLimit(G4double cut)
{
  if (verboseLevel > 0)
    {
      G4cout << "Tst14PhysicsList - Electron cut in energy = " 
	     << cut * MeV << " MeV" << G4endl;
    }  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(cut,1e5);
}

void Tst14PhysicsList::SetGammaCut(G4double value)
{
  //ResetCuts();
  cutForGamma = value;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}


void Tst14PhysicsList::SetElectronCut(G4double value)
{
 // ResetCuts();
  cutForElectron = value;  
  SetParticleCuts(cutForElectron, G4Electron::Electron());

}


void Tst14PhysicsList::SetCuts()
{
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
}



