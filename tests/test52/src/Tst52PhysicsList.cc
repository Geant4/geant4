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
// $Id: Tst52PhysicsList.cc,v 1.2.2.1 2007-12-10 16:34:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: guatelli@ge.infn.it
// -------------------------------------------------------------------

#include "Tst52PhysicsList.hh"
#include "Tst52PhysicsListMessenger.hh"
#include "Tst52Particles.hh"
#include "Tst52PhotonStandard.hh"
#include "Tst52PhotonEPDL.hh"
#include "Tst52PhotonPenelope.hh"
#include "Tst52ElectronStandard.hh"
#include "Tst52ElectronEEDL.hh"
#include "Tst52ElectronPenelope.hh"
#include "Tst52PositronStandard.hh"
#include "Tst52PositronPenelope.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyIonisation.hh"

Tst52PhysicsList::Tst52PhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false)
                                      
{
  defaultCutValue = 0.001 * mm; 
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  SetVerboseLevel(1);

  messenger = new Tst52PhysicsListMessenger(this);
 
  RegisterPhysics( new Tst52Particles("particles") );
  electron_physics_standard = 0;
  electron_physics_eedl = 0;
  electron_physics_penelope = 0;
  positron_physics_standard = 0;
  positron_physics_penelope = 0;
  photon_physics_epdl =0;
}


Tst52PhysicsList::~Tst52PhysicsList()
{
  delete messenger;
}


void Tst52PhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register standard processes for photons
  if (name == "photon-standard") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst52PhotonStandard(name) );
	  photonIsRegistered = true;
	  
	}
    }
  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
          photon_physics_epdl =  new Tst52PhotonEPDL(name);
	  RegisterPhysics(photon_physics_epdl );
	  photonIsRegistered = true;
       
	}
    } 
  // Register processes a' la Penelope for photons
  if (name == "photon-penelope")
    {
      if (photonIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new Tst52PhotonPenelope(name) );
	  photonIsRegistered = true; 
	}
    }

  // Register standard processes for electrons
  if (name == "electron-standard") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  electron_physics_standard =  new Tst52ElectronStandard(name);

	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics(electron_physics_standard);
	  electronIsRegistered = true; 
	  electron_value = name;
          
	}
    }

  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  electron_physics_eedl =  new Tst52ElectronEEDL(name);
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics(electron_physics_eedl);
	  electronIsRegistered = true;
	  electron_value = name;
          
	}
    } 

  // Register processes a' la Penelope for electrons
  if (name == "electron-penelope")
    {
      if (electronIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  electron_physics_penelope =  new Tst52ElectronPenelope(name);
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics(electron_physics_penelope);
	  electronIsRegistered = true;
	  electron_value = name;
	}
    }
  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;

	  positron_physics_standard =  new Tst52PositronStandard(name);
	  RegisterPhysics(positron_physics_standard);
	  positronIsRegistered = true;
	  positron_value = name;
	}
    }


  if (name == "positron-penelope")
    {
      if (positronIsRegistered) 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name 
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "Tst52PhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  positron_physics_penelope =  new Tst52PositronPenelope(name);
	  RegisterPhysics(positron_physics_penelope);
	  positronIsRegistered = true;
	  positron_value = name;
	  
	}
    }


  if (electronIsRegistered && positronIsRegistered && photonIsRegistered)
    {
      G4cout << "PhysicsList for electron, positron and photon registered" << G4endl;
    }
}

void Tst52PhysicsList::SetParticleCut(G4double value)
{ 
  ResetCuts();
  cutForGamma = value;
  cutForElectron = value;
  // defaultCutValue = value; 
}

void Tst52PhysicsList::SetCuts()
{ 
  G4double lowlimit = 250.*eV;
  G4ProductionCutsTable::GetProductionCutsTable()
    ->SetEnergyRange(lowlimit, 100.*GeV);

  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
  // G4VUserPhysicsList::SetCutsWithDefault();
  DumpCutValuesTable();
}

void Tst52PhysicsList::SetFacRange(G4double value)
{ 
  if (electronIsRegistered == true)
    {
      if (electron_value == "electron-standard")  
	electron_physics_standard ->SetFacRange(value);
     
      if (electron_value == "electron-eedl") 
	electron_physics_eedl -> SetFacRange(value);

      if (electron_value == "electron-penelope") 
	electron_physics_penelope -> SetFacRange(value);
    }
  else G4cout <<"Activate the electron physics processes before!!!" << G4cout;

  if (positronIsRegistered == true)
    {
      if (positron_value == "positron-standard")  
	positron_physics_standard -> SetFacRange(value);
     
      if (positron_value == "positron-penelope")  
	positron_physics_penelope -> SetFacRange(value);  
    } else G4cout <<"Activate positron physics processes before!!!" << G4cout;
}




