// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserPhysicsListMessenger.cc,v 1.1 1999-01-07 16:14:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//
//  G4UserPhysicsListMessenger.cc
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
// ------------------------------------------------------------

#include "G4UserPhysicsListMessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

G4UserPhysicsListMessenger::G4UserPhysicsListMessenger(G4VUserPhysicsList* pParticleList):theParticleList(pParticleList)
{
  // /run/particle    directory
  theDirectory = new G4UIdirectory("/run/particle/");
  theDirectory->SetGuidance("Commands for G4VUserPhysicsList.");

  // /run/particle/Verbose command
  verboseCmd = new G4UIcmdWithAnInteger("/run/particle/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level of G4VUserPhysicsList.");
  verboseCmd->SetGuidance(" 0 : Silent (default)");
  verboseCmd->SetGuidance(" 1 : Display warning messages");
  verboseCmd->SetGuidance(" 2 : Display more");
  verboseCmd->SetParameterName("level",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("level >=0 && level <=2");

  // /run/particle/setCut command
  setCutCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setCut",this);
  setCutCmd->SetGuidance("Set default cut value ");
  setCutCmd->SetParameterName("cut",false);
  setCutCmd->SetDefaultValue(1.0);
  setCutCmd->SetRange("cut >0.0");
  setCutCmd->SetUnitCandidates("m cm mm microm");
  setCutCmd->SetDefaultUnit("mm");
  setCutCmd->AvailableForStates(PreInit,Idle);

  // /run/particle/DumpList command
  dumpListCmd = new G4UIcmdWithoutParameter("/run/particle/dumpList",this);
  dumpListCmd->SetGuidance("Dump List of particles in G4VUserPhysicsList. ");


  // /run/particle/DumpCutValues command
  dumpCutValuesCmd = new G4UIcmdWithAString("/run/particle/dumpCutValues",this);
  dumpCutValuesCmd->SetGuidance("Dump cut value information  ");
  dumpCutValuesCmd->SetGuidance("Enter particle name ");
  dumpCutValuesCmd->SetGuidance("    enter all for all particles");
  dumpCutValuesCmd->SetParameterName("particle", true);
  dumpCutValuesCmd->SetDefaultValue("table");
  dumpCutValuesCmd->AvailableForStates(Idle,GeomClosed,EventProc);

  // /run/particle/addProcManager command
  addProcManCmd = new G4UIcmdWithAString("/run/particle/addProcManager", this);
  addProcManCmd->SetGuidance("add process manager");
  addProcManCmd->SetParameterName("particleType", true);
  addProcManCmd->SetDefaultValue("");
  addProcManCmd->AvailableForStates(Init,Idle,GeomClosed,EventProc);
}

G4UserPhysicsListMessenger::~G4UserPhysicsListMessenger()
{
  delete setCutCmd; 
  delete verboseCmd;
  delete dumpListCmd;
  delete dumpCutValuesCmd;
  delete addProcManCmd;
  delete theDirectory;
}

void G4UserPhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==setCutCmd ){
    G4double newCut = setCutCmd->GetNewDoubleValue(newValue); 
    theParticleList->SetDefaultCutValue(newCut);
  } else if( command==verboseCmd ) {
    theParticleList->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)); 

  } else if( command==dumpListCmd ){
    theParticleList->DumpList();
  } else if( command==dumpCutValuesCmd ){ 
    if (newValue == "table") {
      theParticleList->DumpCutValuesTable();
    } else {
      theParticleList->DumpCutValues(newValue);
    }
  }  else if( command == addProcManCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == NULL) return;
    if (particle->GetProcessManager() != NULL) return;
    theParticleList->AddProcessManager(particle);
  }
}

G4String G4UserPhysicsListMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==setCutCmd ){
    cv = setCutCmd->ConvertToString( theParticleList->GetDefaultCutValue(), "mm" );
    
  } else if( command==verboseCmd ){
    cv = verboseCmd->ConvertToString(theParticleList->GetVerboseLevel());

  }  else if( command== addProcManCmd ){
    // set candidate list
    G4String candidates("none");
    G4ParticleTable::G4PTblDicIterator *piter = (G4ParticleTable::GetParticleTable())->GetIterator();
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    addProcManCmd->SetCandidates(candidates);   
    cv = "";
  }
   
  return cv;
}








