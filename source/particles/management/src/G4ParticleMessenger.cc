// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleMessenger.cc,v 1.2 1999-04-13 08:00:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <iomanip.h>                  // Include from 'system'

#include "G4ParticleMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ParticleTable.hh"
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
  thisDirectory->SetGuidance("Paricle control commands.");

  //Commnad   /particle/select
  selectCmd = new G4UIcmdWithAString("/particle/select",this);
  selectCmd->SetGuidance("Select particle ");
  selectCmd->SetDefaultValue("none");
  selectCmd->SetParameterName("particle name", false);

  //Commnad   /particle/list
  listCmd = new G4UIcmdWithAString("/particle/list",this);
  listCmd->SetGuidance("List name of particles.");
  listCmd->SetGuidance(" all(default)/lepton/baryon/meson");
  listCmd->SetParameterName("particle type", true);
  listCmd->SetDefaultValue("all");
  listCmd->SetCandidates("all lepton baryon meson");

  //Commnad   /particle/find 
  findCmd = new G4UIcmdWithAnInteger("/particle/find",this);
  findCmd->SetGuidance("Find particle by encoding");
  findCmd->SetDefaultValue(0);
  findCmd->SetParameterName("encoding", false);

  currentParticle = 0;

  //UI messenger for Particle Properties
  fParticlePropertyMessenger = new G4ParticlePropertyMessenger(theParticleTable);

}

G4ParticleMessenger::~G4ParticleMessenger()
{
  delete fParticlePropertyMessenger;
  delete listCmd; 
  delete selectCmd;
  delete thisDirectory;
  delete findCmd;
}


void G4ParticleMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command==listCmd ){
    //Commnad   /particle/List
    G4int counter = 0;
    G4ParticleTable::G4PTblDicIterator *piter = theParticleTable->GetIterator();
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      newValues.toLower();
      if ((newValues=="all") || (newValues==particle->GetParticleType())) {
        G4cout << setw(19) << particle->GetParticleName();
	if ((counter++)%4 == 3) {
	  G4cout << endl;
	} else {
	  G4cout << ",";
	}
      }
    }
    G4cout << endl;
    if (counter == 0) G4cout << newValues << " is not found " << endl;

  } else if( command==selectCmd ){
    //Commnad   /particle/select
    currentParticle = theParticleTable->FindParticle(newValues);
    if(currentParticle == 0) {
      G4cout << "Unknown particle [" << newValues << "]. Command ignored." << endl;
    }   
  } else if( command==findCmd ){
    //Commnad   /particle/find
    G4ParticleDefinition* tmp = theParticleTable->FindParticle( findCmd->GetNewIntValue(newValues));
    if(tmp == 0) {
      G4cout << "Unknown particle [" << newValues << "]. Command ignored." << endl;
    } else {
      G4cout << tmp->GetParticleName() << endl;
      tmp->DumpTable();
    }
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
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    selectCmd->SetCandidates((const char *)(candidates));   

    static G4String noName("none");
    // current value
    if(currentParticle == 0) {
      // no particle is selected. return null 
      return noName;
    } else {
      return currentParticle->GetParticleName();
    }
  }
  return "";
}






