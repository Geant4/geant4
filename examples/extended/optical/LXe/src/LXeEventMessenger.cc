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





