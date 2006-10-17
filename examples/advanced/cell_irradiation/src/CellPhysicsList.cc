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
//    **************************************
//    *                                    *
//    *        CellPhysicsList.cc          *
//    *                                    *
//    **************************************
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Re-designed for modular Physics List
//
// -------------------------------------------------------------------

#include "CellPhysicsList.hh"
#include "CellPhysicsListMessenger.hh"
#include "CellParticles.hh"
#include "CellPhotonEPDL.hh"
#include "CellElectronEEDL.hh"
#include "CellPositronStandard.hh"
#include "CellProtonICRU49.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

CellPhysicsList::CellPhysicsList(): G4VModularPhysicsList(),
				      electronIsRegistered(false), 
				      positronIsRegistered(false),
				      photonIsRegistered(false), 
                                      protonIsRegistered(false)
{
  // Cut = threshold of pruduction of secondary particles
  // Secondary particles are generated if their range would be bigger
  // than the cut; otherwise the competent energy is considered
  // deposited locally

  defaultCutValue = 0.1 * mm;
  SetVerboseLevel(1);

  // The messenger allows to change physics parameters interactively
  messenger = new CellPhysicsListMessenger(this);
 
  // The particles of the experimental set-up are instantiated
  RegisterPhysics( new CellParticles("particles") );

  // Activate the physics process for every particle 
  AddPhysicsList("photon-epdl");
  AddPhysicsList("electron-eedl");
  AddPhysicsList("positron-standard");
  AddPhysicsList("proton-ICRU49");
}

CellPhysicsList::~CellPhysicsList()
{
  delete messenger;
}

void CellPhysicsList::AddPhysicsList(const G4String& name)
{

  G4cout << "Adding PhysicsList chunk " << name << G4endl;

  // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new CellPhotonEPDL(name) );
	  photonIsRegistered = true;
	}
   } 
 
  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new CellElectronEEDL(name) );
	  electronIsRegistered = true;
	}
   } 

  // Register standard processes for positrons
  if (name == "positron-standard") 
    {
      if (positronIsRegistered) 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- positron List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new CellPositronStandard(name) );
	  positronIsRegistered = true;
	}
    }

// Proton - Low Energy, ICRU49 parameterization 

 if (name == "proton-ICRU49") 
    {
      if (protonIsRegistered) 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- proton e.m. List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "CellPhysicsList::AddPhysicsList: " << name << " is registered" << G4endl;
	  RegisterPhysics( new CellProtonICRU49(name) );
	  protonIsRegistered = true;
	}
    }


  if (electronIsRegistered && 
      positronIsRegistered && 
      photonIsRegistered &&
      protonIsRegistered)
    {
      G4cout << "EM PhysicsList for electron, positron, photon, proton and ions registered" << G4endl;
    }
}

void CellPhysicsList::SetParticleCut(G4double value)
{
  defaultCutValue = value; 
}

void CellPhysicsList::SetCuts()
{
  // This method setd the default CutValue as threshold
  // of production of secondary particles
 G4VUserPhysicsList::SetCutsWithDefault();

 if (verboseLevel>0) DumpCutValuesTable();
}










