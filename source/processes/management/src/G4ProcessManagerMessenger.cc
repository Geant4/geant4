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
// $Id: G4ProcessManagerMessenger.cc 108536 2018-02-16 09:20:49Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ProcessManagerMessenger.cc
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ProcessManagerand UI.
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov. 1997  H. Kurashige   : fixed bugs 
//    08 jan. 1998  H. Kurashige   : new UIcommnds 
//    02 June 2006  M. Maire       : add physicsModified in activate/inactivate
//
//---------------------------------------------------------------


#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManagerMessenger.hh"
#include "G4ios.hh"                 // Include from 'system'
#include <iomanip>                  // Include from 'system'

#include <sstream>

G4ProcessManagerMessenger::G4ProcessManagerMessenger(G4ParticleTable* pTable)
                        :theParticleTable(pTable),
			 currentParticle(0),
			 currentProcess(0),
			 theManager(0),
                         theProcessList(0)
{
  if ( theParticleTable == 0) theParticleTable = G4ParticleTable::GetParticleTable();

  //Commnad   /particle/process
  thisDirectory = new G4UIdirectory("/particle/process/");
  thisDirectory->SetGuidance("Process Manager control commands.");

  //Commnad   /particle/process/dump
  dumpCmd = new G4UIcmdWithAnInteger("/particle/process/dump",this);
  dumpCmd->SetGuidance("dump process manager or process information");
  dumpCmd->SetGuidance("  dump [process index]");
  dumpCmd->SetGuidance("   process index: -1 for process manager");
  dumpCmd->SetParameterName("index", true);
  dumpCmd->SetDefaultValue(-1);

  //Commnad   /particle/process/verbose
  verboseCmd = new G4UIcommand("/particle/process/verbose",this);
  verboseCmd->SetGuidance("Set Verbose Level for Process or Process Manager");
  verboseCmd->SetGuidance("  Verbose [Verbose] [process index]");
  verboseCmd->SetGuidance("   process index: -1 for process manager");
  G4UIparameter* param = new G4UIparameter("Verbose",'i',true);
  param->SetDefaultValue(1);
  verboseCmd->SetParameter(param);
  param = new G4UIparameter("index",'i',true);
  param->SetDefaultValue(-1);
  verboseCmd->SetParameter(param);
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  //Commnad   /particle/process/activate
  activateCmd = new G4UIcmdWithAnInteger("/particle/process/activate",this);
  activateCmd->SetGuidance("Activate process  ");
  activateCmd->SetGuidance(" Activate [process index]");
  activateCmd->SetParameterName("index", false);
  activateCmd->SetDefaultValue(0);
  activateCmd->SetRange("index >=0");
  activateCmd->AvailableForStates(G4State_Idle);

  //Commnad   /particle/process/inactivate
  inactivateCmd = new G4UIcmdWithAnInteger("/particle/process/inactivate",this);
  inactivateCmd->SetGuidance("Inactivate process  ");
  inactivateCmd->SetGuidance(" inactivate [process index]");
  inactivateCmd->SetParameterName("index", false);
  inactivateCmd->SetDefaultValue(0);
  inactivateCmd->SetRange("index >=0");
  inactivateCmd->AvailableForStates(G4State_Idle);

}

G4ProcessManagerMessenger::~G4ProcessManagerMessenger()
{
  delete activateCmd; 
  delete inactivateCmd; 
  delete verboseCmd;
  delete dumpCmd;
  delete thisDirectory;
}

G4ParticleDefinition* G4ProcessManagerMessenger::SetCurrentParticle()
{
  // set currentParticle pointer
  // get particle name by asking G4ParticleMessenger
  G4String particleName = G4UImanager::GetUIpointer()->GetCurrentStringValue("/particle/select");

  currentParticle = theParticleTable->FindParticle(particleName);
  if (currentParticle == 0) {
    theManager = 0;
    G4cout << "G4ProcessManagerMessenger::SetCurrentParticle() ";
    G4cout << particleName << " not found " << G4endl;
  } else {
    theManager = currentParticle->GetProcessManager();
    theProcessList = theManager->GetProcessList();
  }
  return currentParticle;
}

void G4ProcessManagerMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4ExceptionDescription ed;
  if (SetCurrentParticle()==0) {
      ed << "Particle is not selected yet !! Command ignored.";
      command->CommandFailed(ed);
      return;
  }
  if( command == dumpCmd ){
    //Commnad   /particle/process/dump
    G4int index = dumpCmd->GetNewIntValue(newValue);
    if (index <0) {
       theManager->DumpInfo();
    } else if ( index < theManager->GetProcessListLength()){
      currentProcess =  (*theProcessList)(index);
      if (currentProcess == 0) {
	ed << " no process at index of " << index
	   << " in the Process Vector";
        command->CommandFailed(ed);
      } else {
	currentProcess->DumpInfo();
      }
    } else {
      ed << " illegal index !!! ";
      command->CommandFailed(ed);
      currentProcess = 0;
    } 
 
  } else if( command==activateCmd ) {
    //Commnad   /particle/process/activate
    theManager->SetProcessActivation(activateCmd->GetNewIntValue(newValue), true);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    
  } else if( command==inactivateCmd ) {
    //Commnad   /particle/process/inactivate
    theManager->SetProcessActivation(inactivateCmd->GetNewIntValue(newValue), false);
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
    
  } else if( command==verboseCmd ) {
    //Commnad   /particle/process/Verbose
    //  inputstream for newValues 
    const char* temp = (const char*)(newValue);
    std::istringstream is((char*)temp);
    G4int Verbose, index;
    is  >>Verbose >>index;
    if (index <0) {
      theManager->SetVerboseLevel(Verbose);
      
    } else if ( index < theManager->GetProcessListLength()){
      currentProcess =  (*theProcessList)(index);
      if (currentProcess == 0) {
	ed << " no process at index of " << index
	   << " in the Process Vector";
        command->CommandFailed(ed);
      } else {
	currentProcess->SetVerboseLevel(Verbose);
      }
    } else {
      ed << " illegal index !!! ";
      command->CommandFailed(ed);
      currentProcess = 0;
    } 
  }
}


G4String G4ProcessManagerMessenger::GetCurrentValue(G4UIcommand * command)
{
  if(SetCurrentParticle() == 0) return "";

  if( command==verboseCmd ){
    //Commnad   /particle/process/Verbose
    return verboseCmd->ConvertToString(theManager->GetVerboseLevel());
  } else {
    return "";
  }   
}





