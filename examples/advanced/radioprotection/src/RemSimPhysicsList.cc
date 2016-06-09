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
// $Id: RemSimPhysicsList.cc,v 1.6 2004/05/27 13:32:46 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Author: Susanna Guatelli

#include "RemSimPhysicsList.hh"
#include "RemSimPhysicsListMessenger.hh"
#include "RemSimParticles.hh"
#include "RemSimPhotonEPDL.hh"
#include "RemSimElectronEEDL.hh"
#include "RemSimPositronStandard.hh"
#include "RemSimIonICRU.hh"
#include "RemSimIonStandard.hh"
#include "RemSimMuonStandard.hh"
#include "RemSimDecay.hh"
#include "RemSimHadronicPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

RemSimPhysicsList::RemSimPhysicsList(): G4VModularPhysicsList(),
					electronIsRegistered(false), 
					positronIsRegistered(false),
					photonIsRegistered(false), 
                                        ionIsRegistered(false),
                                        hadronicIsRegistered(false),
                                        decayIsRegistered(false),
                                        muonIsRegistered(false)
{
  defaultCutValue = 1.* mm;
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

   // Register LowE-EPDL processes for photons
  if (name == "photon-epdl") 
    {
      if (photonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- photon List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimPhotonEPDL(name) );
	  photonIsRegistered = true;
	}
    } 

  // Register LowE-EEDL processes for electrons
  if (name == "electron-eedl") 
    {
      if (electronIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- electron List already existing"                  << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
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
		 << " cannot be registered ---- positron List already existing"                  << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimPositronStandard(name) );
	  positronIsRegistered = true;
	}
    }

  // Register LowEnergy processes - ICRU Model - for p, alpha, ions. 
  if (name == "ion-ICRU") 
    {
      if (ionIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ----ion e.m. List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimIonICRU(name) );
	  ionIsRegistered = true;
	}
    }

 if (name == "ion-standard") 
    {
      if (ionIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ----ion e.m. List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimIonStandard(name) );
	  ionIsRegistered = true;
	}
    }


  // Register hadronic process for p and alpha particle
  if (name == "hadronic") 
    {
      if (hadronicIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- hadronic physics List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimHadronicPhysics(name) );
	  hadronicIsRegistered = true;
	}
    }

if (name == "muon-standard") 
    {
      if (muonIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered -- muon e.m. List already existing" 
                 << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimMuonStandard(name) );
	  muonIsRegistered = true;
	}
    }

if (name == "decay") 
    {
    if (decayIsRegistered) 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name  
		 << " cannot be registered ---- decay  physics  List already existing" << G4endl;
	} 
      else 
	{
	  G4cout << "RemSimPhysicsList::AddPhysicsList: " << name 
                 << " is registered" << G4endl;
	  RegisterPhysics( new RemSimDecay(name) );
	  decayIsRegistered = true;
	}
    }

  if (electronIsRegistered && positronIsRegistered && photonIsRegistered && 
      ionIsRegistered)
    {
      G4cout << "The electromagnetic processes are registered" << G4endl;
    
  
  if (muonIsRegistered && decayIsRegistered && hadronicIsRegistered) 
    {
      G4cout << "The hadronic physics is registered for p and alpha (up to 100 GeV), the e.m. physics for muons is registered, the decay is registered" 
             << G4endl;
    }
    }
}

void RemSimPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();
}










