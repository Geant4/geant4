// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessManagerMessenger.cc,v 1.2 1999-04-13 09:48:01 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <iomanip.h>                  // Include from 'system'

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

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

  //Commnad   /particle/process/Verbose
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

  //Commnad   /particle/process/Activate
  activateCmd = new G4UIcmdWithAnInteger("/particle/process/activate",this);
  activateCmd->SetGuidance("Activate process  ");
  activateCmd->SetGuidance(" Activate [process index]");
  activateCmd->SetParameterName("index", false);
  activateCmd->SetDefaultValue(0);
  activateCmd->SetRange("index >=0");

  //Commnad   /particle/process/inactivate
  inactivateCmd = new G4UIcmdWithAnInteger("/particle/process/inactivate",this);
  inactivateCmd->SetGuidance("Inactivate process  ");
  inactivateCmd->SetGuidance(" inactivate [process index]");
  inactivateCmd->SetParameterName("index", false);
  inactivateCmd->SetDefaultValue(0);
  inactivateCmd->SetRange("index >=0");
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
    G4cout << particleName << " not found " << endl;
  } else {
    theManager = currentParticle->GetProcessManager();
    theProcessList = theManager->GetProcessList();
  }
  return currentParticle;
}

void G4ProcessManagerMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if (SetCurrentParticle()==0) {
      G4cout << "Particle is not selected yet !! Command ignored." << endl;
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
	G4cout << " no process at index of " << index;
	G4cout << "in the Process Vector" << endl;
      } else {
	currentProcess->DumpInfo();
      }
    } else {
      G4cout << " illegal index !!! " << endl;
      currentProcess = 0;
    } 
 
  } else if( command==activateCmd ) {
    //Commnad   /particle/process/Activate
    theManager->SetProcessActivation(activateCmd->GetNewIntValue(newValue), true);

  } else if( command==inactivateCmd ) {
    //Commnad   /particle/process/inactivate
    theManager->SetProcessActivation(inactivateCmd->GetNewIntValue(newValue), false);
  } else if( command==verboseCmd ) {
    //Commnad   /particle/process/Verbose
    //  inputstream for newValues 
    const char* temp = (const char*)(newValue);
    istrstream is((char*)temp);
    G4int Verbose, index;
    is  >>Verbose >>index;
    if (index <0) {
      theManager->SetVerboseLevel(Verbose);
      
    } else if ( index < theManager->GetProcessListLength()){
      currentProcess =  (*theProcessList)(index);
      if (currentProcess == 0) {
	G4cout << " no process at index of " << index;
	G4cout << "in the Process Vector" << endl;
      } else {
	currentProcess->SetVerboseLevel(Verbose);
      }
    } else {
      G4cout << " illegal index !!! " << endl;
      currentProcess = 0;
    } 
  }
}


G4String G4ProcessManagerMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String returnValue('\0');
  if(SetCurrentParticle() == 0) {
    // no particle is selected. return null strings
    return returnValue;
  }

  char line[255];
  ostrstream os(line,255);
  
  if( command==verboseCmd ){
    //Commnad   /particle/process/Verbose
    os << theManager->GetVerboseLevel() << '\0';
    returnValue = G4String(line);
  } 
  return returnValue;
}





