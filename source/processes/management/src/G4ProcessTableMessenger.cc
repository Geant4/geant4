// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessTableMessenger.cc,v 1.2 1999-04-13 09:48:04 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4ProcessTableMessenger.cc
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ProcessTable and UI.
//
//
//  History:
//    15 Aug. 1998, H. Kurashige  
//---------------------------------------------------------------

#include "G4ProcessTableMessenger.hh"

#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"                 
#include <iomanip.h>               
#include <rw/ctoken.h>               
#include <rw/rstream.h>               

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

/////////////////////////////////////////
G4int G4ProcessTableMessenger::NumberOfProcessType = 10;

//////////////////////////
G4ProcessTableMessenger::G4ProcessTableMessenger(G4ProcessTable* pTable)
                        :theProcessTable(pTable), 
			 currentProcessTypeName("all"),
			 currentProcessName("all"),
			 currentParticleName("all")
{ 
  //Commnad   /particle/process
  thisDirectory = new G4UIdirectory("/process/");
  thisDirectory->SetGuidance("Process Table control commands.");


  //Commnad   /particle/process/list
  listCmd = new G4UIcmdWithAString("/process/list",this);
  listCmd->SetGuidance("List up process names");
  listCmd->SetGuidance("  list [type] ");
  listCmd->SetGuidance("    type: process type [all:for all proceeses]");
  listCmd->SetParameterName("type", true);
  listCmd->SetDefaultValue("all");
  SetNumberOfProcessType();
 
  G4String candidates("all");
  for (G4int idx = 0; idx < NumberOfProcessType ; idx ++ ) {
    candidates += " " + 
      G4VProcess::GetProcessTypeName(G4ProcessType(idx));
  }
  listCmd->SetCandidates((const char*)(candidates));

  //Commnad   /particle/process/Verbose
  verboseCmd = new G4UIcmdWithAnInteger("/process/verbose",this);
  verboseCmd->SetGuidance("Set Verbose Level for Process Table");
  verboseCmd->SetGuidance("  verbose [level]");
  verboseCmd->SetGuidance("   level: verbose level");
  verboseCmd->SetParameterName("verbose", true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose >=0");

  //Commnad   /particle/process/dump
  dumpCmd = new G4UIcommand("/process/dump",this);
  dumpCmd->SetGuidance("Dump process information");
  dumpCmd->SetGuidance(" dump name [particle]");
  dumpCmd->SetGuidance("   name:     process name or type name");
  dumpCmd->SetGuidance("   particle: particle name [all: for all particles]");
  G4UIparameter* param = new G4UIparameter("procName",'s',false);
  dumpCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',true);
  param->SetDefaultValue("all");
  dumpCmd->SetParameter(param);

  //Commnad   /particle/process/activate
  activateCmd = new G4UIcommand("/process/activate",this);
  activateCmd->SetGuidance("Activate processes  ");
  activateCmd->SetGuidance(" Activate  name [particle]");
  activateCmd->SetGuidance("   name:     process name or type name");
  activateCmd->SetGuidance("   particle: particle name [all: for all particles]");
  param = new G4UIparameter("procName",'s',false);
  activateCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',true);
  param->SetDefaultValue("all");
  activateCmd->SetParameter(param);
  
  //Commnad   /particle/process/inactivate
  inactivateCmd = new G4UIcommand("/process/inactivate",this);
  inactivateCmd->SetGuidance("Inactivate process  ");
  inactivateCmd->SetGuidance("Inactivate processes  ");
  inactivateCmd->SetGuidance(" Inactivate  name [particle]");
  inactivateCmd->SetGuidance("   name:     process name or type name");
  inactivateCmd->SetGuidance("   particle: particle name [all: for all particles]");
  param = new G4UIparameter("procName",'s',false);
  inactivateCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',true);
  param->SetDefaultValue("all");
  inactivateCmd->SetParameter(param);
}

//////////////////
G4ProcessTableMessenger::~G4ProcessTableMessenger()
{
  delete activateCmd; 
  delete inactivateCmd; 
  delete verboseCmd;
  delete dumpCmd;
  delete listCmd;
  delete thisDirectory;
}

///////////////
void G4ProcessTableMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4ProcessTable::G4ProcNameVector* procNameVector 
                         = theProcessTable->GetNameList(); 
  G4ProcessTable::G4ProcTableVector* pTblVector 
                         = theProcessTable->GetProcTableVector();
  G4int idx;

  G4ProcessVector* tmpVector;
  G4int type = -1;

  if( command == listCmd ){
    //Commnad  /process/list
    type = -1;
    if (newValue == "all") {	
      currentProcessTypeName = newValue;
    } else {
      type  = GetProcessType(newValue);
      if (type <0) {
	G4cout << " illegal type !!! " << endl;
      } else {
	currentProcessTypeName = newValue;
      }
    }    
    G4int counter = 0;
    for (idx = 0; idx < procNameVector->length(); idx++) {
      G4String name = (*procNameVector)(idx);
      tmpVector = theProcessTable->FindProcesses(name);
      if ( (type <0) || ( ((*tmpVector)(0)->GetProcessType()) == type) ) {
        if ( counter%4 != 0) G4cout << ",";
	G4cout << setw(19) << name;
	if ((counter++)%4 == 3) {
          G4cout << endl;
        }
      }
    }
    G4cout << endl;
    delete tmpVector;
    //Commnad  /process/list

  } else if( command==verboseCmd ) {
    //Commnad   /process/verbose
     theProcessTable->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));
    //Commnad   /process/verbose

  } else {
    RWCTokenizer next( newValue );

    // check 1st argument
    currentProcessName = G4String(next());
    G4bool isNameFound = false;
    G4bool isProcName = false; 
    if ( procNameVector->contains(currentProcessName) ){
      isNameFound = true;
      isProcName  = true; 
    } else {
      type  = GetProcessType(currentProcessName);
      if (type >=0) {
	isNameFound = true;
      } else {
	// no processes with specifed name
	G4cout << " illegal process (or type) name " << endl;
	currentProcessName = "";
	return;
      }
    }
  
    // check 2nd argument
    currentParticleName = G4String(next());
    G4bool isParticleFound = false;
    G4ParticleDefinition* currentParticle = 0;
    if ( currentParticleName == "all" ) {
      isParticleFound = true;

    } else {
      isParticleFound = G4ParticleTable::GetParticleTable()->contains(currentParticleName);
      if (isParticleFound) {
	currentParticle = G4ParticleTable::GetParticleTable()->FindParticle(currentParticleName);
      }

    }

    if ( !isParticleFound ) {
      // no particle with specifed name
      G4cout << " illegal particle name " << endl;
      currentParticleName = "";
      return;
    }
        
    if( command==dumpCmd ) {
      // process/dump
      if (isProcName) {
	tmpVector = theProcessTable->FindProcesses(currentProcessName);
      } else {
	tmpVector = theProcessTable->FindProcesses(G4ProcessType(type));
      }
      for (G4int i=0; i<tmpVector->length(); i++) {
	theProcessTable->DumpInfo( (*tmpVector)(i), currentParticle );
      }
      delete tmpVector;
      // process/dump

    } else if ( (command==activateCmd) || (command==inactivateCmd)) {
      // process/activate , inactivate
      G4bool fActive = (command==activateCmd);
      if (isProcName) {
	if ( currentParticle == 0 ) {
	  theProcessTable->SetProcessActivation(currentProcessName, 
						fActive);
	} else {
	  theProcessTable->SetProcessActivation(currentProcessName,
						currentParticle,
						fActive);
	}
      } else {
	if ( currentParticle == 0 ) {
	  theProcessTable->SetProcessActivation(G4ProcessType(type),
						fActive);
	} else {
	  theProcessTable->SetProcessActivation(G4ProcessType(type),
						currentParticle,
						fActive);
	}
      }
      //  process/activate , inactivate
    } 
  }
}


//////////////////
G4String G4ProcessTableMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4ProcessTable::G4ProcNameVector* procNameVector 
                         = theProcessTable->GetNameList(); 

  G4String candidates;
  G4String returnValue('\0');

  char line[255];
  ostrstream os(line,255);
  G4UIparameter * param; 

  G4int idx; 

  if( command==verboseCmd ){
    //Commnad   /process/verbose
    os << theProcessTable->GetVerboseLevel() << '\0';
    returnValue = G4String(line);

  } else if ( command==listCmd ){
    //Commnad   /process/list
    candidates = "all";
    for (idx = 0; idx < NumberOfProcessType ; idx ++ ) {
      candidates += " " + 
	       G4VProcess::GetProcessTypeName(G4ProcessType(idx));
    }
    listCmd->SetCandidates((const char*)(candidates));
    returnValue =  currentProcessTypeName;

  } else {
    //Commnad   /process/dump, activate, inactivate
    // process name 
    param = command->GetParameter(0);
    candidates = "";
    for (idx = 0; idx < procNameVector->length() ; idx ++ ) {
      candidates += " " + (*procNameVector)(idx);
    }
    param->SetParameterCandidates((const char*)(candidates));
    // particle name
    param = command->GetParameter(1);
    candidates = "all";
    G4ParticleTable::G4PTblDicIterator *piter 
                        = G4ParticleTable::GetParticleTable()->GetIterator();
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    param->SetParameterCandidates((const char*)(candidates));

    returnValue =  currentProcessName + " " + currentParticleName;

  }

  return returnValue;
}

/////////////////
G4String G4ProcessTableMessenger::GetProcessTypeName(G4ProcessType aType) const
{
  return G4VProcess::GetProcessTypeName(aType);
}

/////////////////
G4int G4ProcessTableMessenger::GetProcessType(const G4String& aTypeName) const
{
  G4int type = -1;
  for (G4int idx = 0; idx < NumberOfProcessType ; idx ++ ) {
    if (aTypeName == G4VProcess::GetProcessTypeName(G4ProcessType(idx)) ) {
      type = idx;
      break;
    }
  }
  return type;
}


/////////////////
void G4ProcessTableMessenger::SetNumberOfProcessType()
{
  G4bool isFoundEndMark = false;
  G4int idx;
  for (idx = 0; idx < 1000 ; idx ++ ) {
    G4String typeName = G4VProcess::GetProcessTypeName(G4ProcessType(idx));
    if ( isFoundEndMark = typeName.contains("---")) break;
  }
  if ( isFoundEndMark ) {
    NumberOfProcessType = idx;
  } else {
    G4Exception("G4ProcessTableMessenger::SetNumberOfProcessType(): No End Mark for  G4VProcess::GetProcessTypeName() ");
  } 
}




