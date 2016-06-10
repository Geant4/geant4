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
//
// $Id: G4ParticlePropertyMessenger.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ParticlePropertyMessenger.cc
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ParticleDefinition and UI.
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov. 1997  H. Kurashige   : fixed bugs 
//    08 jan. 1998  H. Kurashige   : new UIcommnds 
//    19 June 1998  H. Kurashige   : adds UnitCategory 
//---------------------------------------------------------------

#include "G4ParticlePropertyMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTableMessenger.hh"
#include "G4ParticlePropertyMessenger.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"                 // Include from 'system'
#include <iomanip>                  // Include from 'system'

G4ParticlePropertyMessenger::G4ParticlePropertyMessenger(G4ParticleTable* pTable)
                        :theParticleTable(pTable),
			 currentParticle(0),
			 fDecayTableMessenger(0)
{
  if ( theParticleTable == 0) theParticleTable = G4ParticleTable::GetParticleTable();
  //Commnad   /particle/property/
  thisDirectory = new G4UIdirectory("/particle/property/");
  thisDirectory->SetGuidance("Paricle Table control commands.");

  //Commnad   /particle/property/dump
  dumpCmd = new G4UIcmdWithoutParameter("/particle/property/dump",this);
  dumpCmd->SetGuidance("dump particle properties.");

  //Command   /particle/property/stable
  stableCmd = new G4UIcmdWithABool("/particle/property/stable",this);
  stableCmd->SetGuidance("Set stable flag.");
  stableCmd->SetGuidance("  false: Unstable   true: Stable");
  stableCmd->SetParameterName("stable",false);
  stableCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);

  //particle/property/lifetime
  lifetimeCmd = new G4UIcmdWithADoubleAndUnit("/particle/property/lifetime",this);
  lifetimeCmd->SetGuidance("Set life time.");
  lifetimeCmd->SetGuidance("Unit of the time can be :");
  lifetimeCmd->SetGuidance(" s, ms, ns (default)");
  lifetimeCmd->SetParameterName("life",false);
  lifetimeCmd->SetDefaultValue(0.0);
  lifetimeCmd->SetRange("life >0.0");
  //lifetimeCmd->SetUnitCategory("Time");
  //lifetimeCmd->SetUnitCandidates("s ms ns");
  lifetimeCmd->SetDefaultUnit("ns");
  lifetimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed);

  // -- particle/property/Verbose ---
  verboseCmd = new G4UIcmdWithAnInteger("/particle/property/verbose",this);
  verboseCmd->SetGuidance("Set Verbose level of particle property.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display warning messages");
  verboseCmd->SetGuidance(" 2 : Display more");
  verboseCmd->SetParameterName("verbose_level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("verbose_level >=0");

  //UI messenger for Decay Table
  fDecayTableMessenger = new G4DecayTableMessenger(theParticleTable);
  
}

G4ParticlePropertyMessenger::~G4ParticlePropertyMessenger()
{
  if (fDecayTableMessenger !=0) delete  fDecayTableMessenger;
  fDecayTableMessenger = 0;

  delete stableCmd; 
  delete verboseCmd;
  delete lifetimeCmd;
  delete dumpCmd;
  delete thisDirectory;
} 

void G4ParticlePropertyMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if (SetCurrentParticle()==0) {
      G4cout << "Particle is not selected yet !! Command ignored." << G4endl;
      return;
  }

  if( command == dumpCmd ){
    //Commnad   /particle/property/dump
    currentParticle->DumpTable();

  } else if (command == lifetimeCmd ) {
    //Commnad   /particle/property/lifetime
    currentParticle->SetPDGLifeTime(lifetimeCmd->GetNewDoubleValue(newValue)); 

  } else if (command == stableCmd ) {
    //Commnad   /particle/property/stable
    if (currentParticle->GetPDGLifeTime()<0.0) {
      G4cout << "Life time is negative! Command ignored." << G4endl; 
    } else if (currentParticle->GetPDGMass()<=0.0) {
      G4cout << "Zero Mass! Command ignored." << G4endl; 
    } else {
      currentParticle->SetPDGStable(stableCmd->GetNewBoolValue(newValue));
    }
  
  } else if( command==verboseCmd ) {
    //Commnad   /particle/property/Verbose
    currentParticle->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); 
  }
}

G4ParticleDefinition* G4ParticlePropertyMessenger::SetCurrentParticle()
{
  // set currentParticle pointer
  
  // get particle name by asking G4ParticleMessenger via UImanager

  G4String particleName = G4UImanager::GetUIpointer()->GetCurrentStringValue("/particle/select");
	
  if (currentParticle != 0 ){
    // check whether selection is changed 
    if (currentParticle->GetParticleName() != particleName) {
      currentParticle = theParticleTable->FindParticle(particleName);
    }
  } else {
    currentParticle = theParticleTable->FindParticle(particleName);
  }
  return currentParticle;
}

G4String G4ParticlePropertyMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String returnValue('\0');

  if (SetCurrentParticle()==0) {
    // no particle is selected. return null 
    return returnValue;
  }

  if( command == stableCmd ){
    //Commnad   /particle/property/stable
    returnValue = stableCmd->ConvertToString( currentParticle->GetPDGStable());

  } else if( command == lifetimeCmd ){
    //Commnad   /particle/property/lifetime
    returnValue = lifetimeCmd->ConvertToString(  currentParticle->GetPDGLifeTime() , "ns" );
    
  } else if( command==verboseCmd ){
   //Commnad   /particle/property/Verbose
     returnValue= verboseCmd->ConvertToString(currentParticle ->GetVerboseLevel());

  } 
  
  return returnValue;
}








