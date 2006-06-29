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
// $Id: RemSimPhysicsList.cc,v 1.10 2006-06-29 16:24:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli

#include "RemSimPhysicsList.hh"
#include "RemSimPhysicsListMessenger.hh"
#include "RemSimParticles.hh"
#include "RemSimPhotonEPDL.hh"
#include "RemSimElectronEEDL.hh"
#include "RemSimPositronStandard.hh"
#include "RemSimIonICRU.hh"
#include "RemSimMuonStandard.hh"
#include "RemSimDecay.hh"
#include "RemSimHadronicBertini.hh"
#include "RemSimHadronicBinary.hh"
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
  //  SetVerboseLevel(1);

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

 //  // Register hadronic process for p and alpha particle, activating the
//   // Binary model

 
 if (name == "hadronic-binary") 
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
	  RegisterPhysics( new RemSimHadronicBinary(name) );
	  hadronicIsRegistered = true;
	}
    }


  // Register hadronic process for p and alpha particle, activating the
  // Bertini model

if (name == "hadronic-bertini") 
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
	  RegisterPhysics( new RemSimHadronicBertini(name) );
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
      G4cout << "The hadronic physics is registered for p (up to 100 GeV) and alpha (up to 10 GeV), the e.m. physics for muons is registered, the decay is registered" 
             << G4endl;
    }
    }
}

void RemSimPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();
}










