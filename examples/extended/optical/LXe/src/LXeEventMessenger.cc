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
#include "LXeEventMessenger.hh"
#include "LXeEventAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeEventMessenger::LXeEventMessenger(LXeEventAction* event)
:LXeEvent(event)
{
  saveThresholdCmd = new G4UIcmdWithAnInteger("/LXe/saveThreshold",this);
  saveThresholdCmd->SetGuidance("Set the photon count threshold for saving the random number seed for an event.");
  saveThresholdCmd->SetParameterName("photons",true);
  saveThresholdCmd->SetDefaultValue(4500);
  saveThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/LXe/eventVerbose",this);
  verboseCmd->SetGuidance("Set the verbosity of event data.");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);

  pmtThresholdCmd = new G4UIcmdWithAnInteger("/LXe/pmtThreshold",this);
  pmtThresholdCmd->SetGuidance("Set the pmtThreshold (in # of photons)");

  forceDrawPhotonsCmd=new G4UIcmdWithABool("/LXe/forceDrawPhotons",this);
  forceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
  forceDrawPhotonsCmd
    ->SetGuidance("(Higher priority than /LXe/forceDrawNoPhotons)");

  forceDrawNoPhotonsCmd=new G4UIcmdWithABool("/LXe/forceDrawNoPhotons",this);
  forceDrawNoPhotonsCmd->SetGuidance("Force no drawing of photons.");
  forceDrawNoPhotonsCmd
    ->SetGuidance("(Lower priority than /LXe/forceDrawPhotons)");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeEventMessenger::~LXeEventMessenger(){
  delete saveThresholdCmd;
  delete verboseCmd;
  delete pmtThresholdCmd;
  delete forceDrawPhotonsCmd;
  delete forceDrawNoPhotonsCmd;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeEventMessenger::SetNewValue(G4UIcommand* command, G4String newValue){ 
  if( command == saveThresholdCmd ){ 
    LXeEvent->SetSaveThreshold(saveThresholdCmd->GetNewIntValue(newValue));
  }
  else if( command == verboseCmd ){
    LXeEvent->SetEventVerbose(verboseCmd->GetNewIntValue(newValue));
  }
  else if( command == pmtThresholdCmd ){
    LXeEvent->SetPMTThreshold(pmtThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == forceDrawPhotonsCmd){
    LXeEvent->SetForceDrawPhotons(forceDrawPhotonsCmd
				  ->GetNewBoolValue(newValue));
  }
  else if(command == forceDrawNoPhotonsCmd){
    LXeEvent->SetForceDrawNoPhotons(forceDrawNoPhotonsCmd
				  ->GetNewBoolValue(newValue));
    G4cout<<"TEST"<<G4endl;
  }
}





