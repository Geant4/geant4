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
// $Id: G4ParticleMessenger.cc 91885 2015-08-10 07:05:56Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ParticleMessenger.cc
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between Particle related classes and UI.
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov. 1997  H. Kurashige   : fixed bugs 
//    08 Jan. 1998  H. Kurashige   : new UIcommnds 
//    08 June 1998, H. Kurashige   : remove fProcessManagerMessenger
//    25 Nov. 1998, H. Kurashige   : add /particle/find
//---------------------------------------------------------------

#include "G4ios.hh"                 // Include from 'system'
#include <iomanip>                  // Include from 'system'

#include "G4ParticleMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticlePropertyMessenger.hh"

G4ParticleMessenger::G4ParticleMessenger(G4ParticleTable* pTable)
{
  // get the pointer to ParticleTable
  if ( pTable == 0) {
    theParticleTable = G4ParticleTable::GetParticleTable();
  } else {
    theParticleTable = pTable;
  }
 
  //Directory   /particle/
  thisDirectory = new G4UIdirectory("/particle/");
  thisDirectory->SetGuidance("Particle control commands.");

  //Commnad   /particle/select
  selectCmd = new G4UIcmdWithAString("/particle/select",this);
  selectCmd->SetGuidance("Select particle ");
  selectCmd->SetDefaultValue("none");
  selectCmd->SetParameterName("particle name", false);
  selectCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //Commnad   /particle/list
  listCmd = new G4UIcmdWithAString("/particle/list",this);
  listCmd->SetGuidance("List name of particles.");
  listCmd->SetGuidance(" all(default)/lepton/baryon/meson/nucleus/quarks");
  listCmd->SetParameterName("particle type", true);
  listCmd->SetDefaultValue("all");
  listCmd->SetCandidates("all lepton baryon meson nucleus quarks");
  listCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //Commnad   /particle/find 
  findCmd = new G4UIcmdWithAnInteger("/particle/find",this);
  findCmd->SetGuidance("Find particle by encoding");
  findCmd->SetDefaultValue(0);
  findCmd->SetParameterName("encoding", false);
  findCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //Commnad   /particle/createAllIon
  createAllIonCmd = new G4UIcmdWithoutParameter("/particle/createAllIon",this);
  createAllIonCmd->SetGuidance("Create All ions (ground state)");
  createAllIonCmd->AvailableForStates(G4State_Idle);
  createAllIonCmd->SetToBeBroadcasted(false);

  //Commnad   /particle/createAllIsomer
  createAllIsomerCmd = new G4UIcmdWithoutParameter("/particle/createAllIsomer",this);
  createAllIsomerCmd->SetGuidance("Create All isomers");
  createAllIsomerCmd->AvailableForStates(G4State_Idle);
  createAllIsomerCmd->SetToBeBroadcasted(false);

  // -- particle/property/Verbose ---
  verboseCmd = new G4UIcmdWithAnInteger("/particle/verbose",this);
  verboseCmd->SetGuidance("Set Verbose level of particle table.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display warning messages");
  verboseCmd->SetGuidance(" 2 : Display more");
  verboseCmd->SetParameterName("verbose_level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("verbose_level >=0");

  currentParticle = 0;

  //UI messenger for Particle Properties
  fParticlePropertyMessenger = new G4ParticlePropertyMessenger(theParticleTable);

}

G4ParticleMessenger::~G4ParticleMessenger()
{
  delete fParticlePropertyMessenger;

  delete listCmd; 
  delete selectCmd;
  delete findCmd;
  delete createAllIonCmd;
  delete createAllIsomerCmd;
  delete verboseCmd;

  delete thisDirectory;
}


void G4ParticleMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command==listCmd ){
    //Commnad   /particle/List
    G4int counter = 0;
    G4ParticleTable::G4PTblDicIterator *piter = theParticleTable->GetIterator();
    piter -> reset();
  
    while( (*piter)() ){// Loop checking, 09.08.2015, K.Kurashige
      G4ParticleDefinition *particle = piter->value();
      if ((newValues=="all") || (newValues==particle->GetParticleType())) {
        G4cout << std::setw(19) << particle->GetParticleName();
	if ((counter++)%4 == 3) {
	  G4cout << G4endl;
	} else {
	  G4cout << ",";
	}
      }
    }
    G4cout << G4endl;
    if (counter == 0) G4cout << newValues << " is not found " << G4endl;
    
   //Command  /particle/select
    // set candidate List
    G4String candidates("none");
    piter -> reset();
    while( (*piter)() ){// Loop checking, 09.08.2015, K.Kurashige
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    selectCmd->SetCandidates((const char *)(candidates));   
 
  } else if( command==selectCmd ){
    //Commnad   /particle/select
    currentParticle = theParticleTable->FindParticle(newValues);
    if(currentParticle == 0) {
      G4cout << "Unknown particle [" << newValues << "]. Command ignored." << G4endl;
    }   

  } else if( command==findCmd ){
    //Commnad   /particle/find
    G4ParticleDefinition* tmp = theParticleTable->FindParticle( findCmd->GetNewIntValue(newValues));
    if(tmp == 0) {
      G4cout << "Unknown particle [" << newValues << "]. Command ignored." << G4endl;
    } else {
      G4cout << tmp->GetParticleName() << G4endl;
      tmp->DumpTable();
    }

  } else if( command==createAllIonCmd ) {
    //Commnad   /particle/createAllIon
    theParticleTable->GetIonTable()->CreateAllIon();

  } else if( command==createAllIsomerCmd ) {
    //Commnad   /particle/createAllIsomer
    theParticleTable->GetIonTable()->CreateAllIsomer();

  } else if( command==verboseCmd ) {
    //Commnad   /particle/verbose
    theParticleTable->SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));
  }    
}

G4String G4ParticleMessenger::GetCurrentValue(G4UIcommand * command)
{
  if( command==selectCmd ){
    //Command  /particle/select
    // set candidate List
    G4String candidates("none");
    G4ParticleTable::G4PTblDicIterator *piter = theParticleTable->GetIterator();
    piter -> reset();
    while( (*piter)() ){// Loop checking, 09.08.2015, K.Kurashige
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    selectCmd->SetCandidates((const char *)(candidates));   

    static const G4String noName("none");
    // current value
    if(currentParticle == 0) {
      // no particle is selected. return null 
      return noName;
    } else {
      return currentParticle->GetParticleName();
    }
  } else if( command==verboseCmd ){
    //Commnad   /particle/verbose
    return verboseCmd->ConvertToString(theParticleTable->GetVerboseLevel());
  }
  return "";
}






