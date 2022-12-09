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
// G4ProcessTableMessenger class implementation
//
// Author: H.Kurashige, 15 August 1998
// --------------------------------------------------------------------

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
#include "G4Tokenizer.hh"
#include <iomanip>               
#include <sstream>

G4ThreadLocal G4int G4ProcessTableMessenger::NumberOfProcessType = 10;

// --------------------------------------------------------------------
G4ProcessTableMessenger::G4ProcessTableMessenger(G4ProcessTable* pTable)
  : theProcessTable(pTable)
{ 
  // Command   /particle/process
  thisDirectory = new G4UIdirectory("/process/");
  thisDirectory->SetGuidance("Process Table control commands.");

  // Command   /particle/process/list
  listCmd = new G4UIcmdWithAString("/process/list",this);
  listCmd->SetGuidance("List up process names");
  listCmd->SetGuidance("  list [type] ");
  listCmd->SetGuidance("    type: process type [all:for all processes]");
  listCmd->SetParameterName("type", true);
  listCmd->SetDefaultValue("all");
  SetNumberOfProcessType();
 
  G4String candidates("all");
  for (G4int idx = 0; idx < NumberOfProcessType ; ++idx )
  {
    candidates += " " + G4VProcess::GetProcessTypeName(G4ProcessType(idx));
  }
  listCmd->SetCandidates((const char*)(candidates));

  // Command   /particle/process/verbose
  verboseCmd = new G4UIcmdWithAnInteger("/process/verbose",this);
  verboseCmd->SetGuidance("Set Verbose Level for Process Table");
  verboseCmd->SetGuidance("  verbose [level]");
  verboseCmd->SetGuidance("   level: verbose level");
  verboseCmd->SetParameterName("verbose", true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose >=0");
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // Command   /particle/process/setVerbose
  procVerboseCmd = new G4UIcommand("/process/setVerbose",this);
  procVerboseCmd->SetGuidance("Set verbose level for processes");
  procVerboseCmd->SetGuidance("  setVerbose level [type or name] ");
  procVerboseCmd->SetGuidance("    level: verbose level ");
  procVerboseCmd->SetGuidance("    name : process name ");
  procVerboseCmd->SetGuidance("    type : process type ");
  procVerboseCmd->SetGuidance("       [all] for all processes ");
  G4UIparameter* param = new G4UIparameter("verbose",'i',false);
  procVerboseCmd->SetParameter(param);
  param = new G4UIparameter("type",'s',true);
  param->SetDefaultValue("all");
  procVerboseCmd->SetParameter(param);
  procVerboseCmd->AvailableForStates(G4State_Idle,G4State_GeomClosed,G4State_EventProc);
 
  // Command   /particle/process/dump
  dumpCmd = new G4UIcommand("/process/dump",this);
  dumpCmd->SetGuidance("Dump process information");
  dumpCmd->SetGuidance(" dump name [particle]");
  dumpCmd->SetGuidance("   name:     process name or type name");
  dumpCmd->SetGuidance("   particle: particle name [all: for all particles]");
  param = new G4UIparameter("procName",'s',false);
  dumpCmd->SetParameter(param);
  param = new G4UIparameter("particle",'s',true);
  param->SetDefaultValue("all");
  dumpCmd->SetParameter(param);
  dumpCmd->AvailableForStates(G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // Command   /process/activate
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
  activateCmd->AvailableForStates(G4State_Idle);
  
  // Command   /process/inactivate
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
  inactivateCmd->AvailableForStates(G4State_Idle);
}

// --------------------------------------------------------------------
G4ProcessTableMessenger::~G4ProcessTableMessenger()
{
  delete activateCmd; 
  delete inactivateCmd; 
  delete verboseCmd;
  delete dumpCmd;
  delete listCmd;
  delete procVerboseCmd;
  delete thisDirectory;
}

// --------------------------------------------------------------------
void
G4ProcessTableMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  G4ProcessTable::G4ProcNameVector* procNameVector 
                         = theProcessTable->GetNameList(); 
  G4int type = -1;
  G4ExceptionDescription ed;

  if( command == listCmd )
  {
    // Command  /process/list
    type = -1;
    if (newValue == "all")
    {
      currentProcessTypeName = newValue;
    }
    else
    {
      type  = GetProcessType(newValue);
      if (type <0)
      {
        G4cout << " illegal type !!! " << G4endl;
      }
      else
      {
        currentProcessTypeName = newValue;
      }
    }    
    G4int counter = 0;
    for (auto itr=procNameVector->cbegin(); itr!=procNameVector->cend(); ++itr)
    {
      G4ProcessVector* tmpVector = theProcessTable->FindProcesses(*itr);
      if ( (type <0) || ( ((*tmpVector)(0)->GetProcessType()) == type) )
      {
        if ( counter%4 != 0) G4cout << ",";
        G4cout << std::setw(19) << *itr;
        if ((counter++)%4 == 3)
        {
          G4cout << G4endl;
        }
      }
      delete tmpVector;
    }
    G4cout << G4endl;
  }
  else if( command==procVerboseCmd )
  {
    // Command  /process/setVerbose
    G4Tokenizer next( newValue );

    // check 1st argument
    G4String tmpS = G4String(next());
    //  inputstream for newValues
    const char* temp = (const char*)(tmpS);
    std::istringstream is((char*)temp);
    G4int level;
    is  >>level;

    // check 2nd argument
    currentProcessTypeName = G4String(next());
    if (currentProcessTypeName.empty()) currentProcessTypeName = "all";
    G4bool isProcName = false;
    G4bool isAll = false;
    type = -1;

    if (currentProcessTypeName == "all")
    {
      isAll = true; 
    }
    else
    {
      type  = GetProcessType(currentProcessTypeName);
      if (type<0)
      {
        isProcName = true;
        currentProcessName = currentProcessTypeName;
        currentProcessTypeName = "";
      }
    }  
    for (auto itr=procNameVector->cbegin(); itr!=procNameVector->cend(); ++itr)
    {
      G4ProcessVector* tmpVector = theProcessTable->FindProcesses(*itr);
      G4VProcess* p = (*tmpVector)(0);
      if ( isAll || 
           (!isProcName && ( p->GetProcessType() == type) ) ||
           ( isProcName && ( p->GetProcessName()== currentProcessName) ) )
      {
        p->SetVerboseLevel(level);
      }
      delete tmpVector;
    }
  }
  else if( command==verboseCmd )
  {
    // Command   /process/verbose
    theProcessTable->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));
  }
  else
  {
    G4Tokenizer next( newValue );

    // check 1st argument
    currentProcessName = G4String(next());
    G4bool isProcName = false; 
    for (auto itr=procNameVector->cbegin(); itr!=procNameVector->cend(); ++itr)
    {
      if ( (*itr) == currentProcessName )
      {
        isProcName = true; 
        break;
      }
    }
    if (!isProcName)
    {
      type  = GetProcessType(currentProcessName);
      if (type <0 )
      {
        // no processes with specified name
        ed << " illegal process (or type) name ["
           << currentProcessName << "]";
        command->CommandFailed(ed);
        currentProcessName = "";
        return;
      }
    }
  
    // check 2nd argument
    currentParticleName = G4String(next());
    G4bool isParticleFound = false;
    G4ParticleDefinition* currentParticle = nullptr;
    if ( currentParticleName == "all" )
    {
      isParticleFound = true;
    }
    else
    {
      isParticleFound = G4ParticleTable::GetParticleTable()
                      ->contains(currentParticleName);
      if (isParticleFound)
      {
        currentParticle = G4ParticleTable::GetParticleTable()
                        ->FindParticle(currentParticleName);
      }
    }

    if ( !isParticleFound )
    {
      // no particle with specified name
      ed << " illegal particle name [" << currentParticleName << "]";
      command->CommandFailed(ed);
      currentParticleName = "";
      return;
    }
        
    if( command==dumpCmd )
    {
      // process/dump
      G4ProcessVector* tmpVector;
      if (isProcName)
      {
        tmpVector = theProcessTable->FindProcesses(currentProcessName);
      }
      else
      {
        tmpVector = theProcessTable->FindProcesses(G4ProcessType(type));
      }
      for (G4int i=0; i<(G4int)tmpVector->length(); ++i)
      {
        theProcessTable->DumpInfo( (*tmpVector)(i), currentParticle );
      }
      delete tmpVector;
    }
    else if ( (command==activateCmd) || (command==inactivateCmd))
    {
      // process/activate , inactivate
      G4bool fActive = (command==activateCmd);
      if (isProcName)
      {
        if ( currentParticle == nullptr )
        {
          theProcessTable->SetProcessActivation(currentProcessName, 
                                                fActive);
        }
        else
        {
          theProcessTable->SetProcessActivation(currentProcessName,
                                                currentParticle,
                                                fActive);
        }
      }
      else
      {
        if ( currentParticle == nullptr )
        {
          theProcessTable->SetProcessActivation(G4ProcessType(type),
                                                fActive);
        }
        else
        {
          theProcessTable->SetProcessActivation(G4ProcessType(type),
                                                currentParticle,
                                                fActive);
        }
      }
      G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    } 
  }
}

// --------------------------------------------------------------------
G4String G4ProcessTableMessenger::GetCurrentValue(G4UIcommand* command)
{
  if( command==verboseCmd )
  {
    // Command   /process/verbose
    return verboseCmd->ConvertToString(theProcessTable->GetVerboseLevel());
  }
  else if ( command==listCmd )
  {
    // Command   /process/list
    return currentProcessTypeName;
  }
  else
  {
    // Command   /process/dump, activate, inactivate
    return   (currentProcessName + " " + currentParticleName);
  }

  return "";
}

// --------------------------------------------------------------------
G4String G4ProcessTableMessenger::GetProcessTypeName(G4ProcessType aType) const
{
  return G4VProcess::GetProcessTypeName(aType);
}

// --------------------------------------------------------------------
G4int G4ProcessTableMessenger::GetProcessType(const G4String& aTypeName) const
{
  G4int type = -1;
  for (G4int idx = 0; idx < NumberOfProcessType ; ++idx )
  {
    if (aTypeName == G4VProcess::GetProcessTypeName(G4ProcessType(idx)) )
    {
      type = idx;
      break;
    }
  }
  return type;
}

// --------------------------------------------------------------------
void G4ProcessTableMessenger::SetNumberOfProcessType()
{
  G4bool isFoundEndMark = false;
  G4int idx;
  for (idx = 0; idx < 1000 ; ++idx )
  {
    G4String typeName = G4VProcess::GetProcessTypeName(G4ProcessType(idx));
    isFoundEndMark = G4StrUtil::contains(typeName, "---");
    if ( isFoundEndMark ) break;
  }
  if ( isFoundEndMark )
  {
    NumberOfProcessType = idx;
  }
  else
  {
    G4Exception("G4ProcessTableMessenger::SetNumberOfProcessType()",
                "ProcMan014", FatalException, "No End Mark");
  } 
}
