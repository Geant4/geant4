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
 
#include "G4HumanPhantomMessenger.hh"
#include "G4HumanPhantomConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "globals.hh"

#include "G4RunManager.hh"

G4HumanPhantomMessenger::G4HumanPhantomMessenger(G4HumanPhantomConstruction* myUsrPhtm)
  :myUserPhantom(myUsrPhtm),bps(false)
{ 
  phantomDir = new G4UIdirectory("/phantom/");
  phantomDir->SetGuidance("Set Your Phantom.");
  
  bpDir = new G4UIdirectory("/bodypart/");
  bpDir->SetGuidance("Add Body Part to Phantom");

  cleanCmd = new G4UIcmdWithoutParameter("/phantom/startNewPhantom",this);
  cleanCmd->SetGuidance("Start a New Phantom and clean previous one.");
  cleanCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  modelCmd = new G4UIcmdWithAString("/phantom/setPhantomModel",this);
  modelCmd->SetGuidance("Set sex of Phantom: MIRD, ORNL or MIX.");
  modelCmd->SetParameterName("phantomModel",true);
  modelCmd->SetDefaultValue("MIRD");
  modelCmd->SetCandidates("MIRD ORNL MIX");
  modelCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  sexCmd = new G4UIcmdWithAString("/phantom/setPhantomSex",this);
  sexCmd->SetGuidance("Set sex of Phantom: Male or Female.");
  sexCmd->SetParameterName("phantomSex",true);
  sexCmd->SetDefaultValue("Female");
  sexCmd->SetCandidates("Male Female");
  sexCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  bodypartCmd = new G4UIcmdWithAString("/bodypart/addBodyPart",this);
  bodypartCmd->SetGuidance("Add a Body Part to Phantom");
  bodypartCmd->SetParameterName("bpName",true);
  bodypartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  endCmd = new G4UIcmdWithoutParameter("/phantom/buildNewPhantom",this);
  endCmd->SetGuidance("Build your Phantom.");
  endCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

G4HumanPhantomMessenger::~G4HumanPhantomMessenger()
{
  delete  cleanCmd;
  delete  modelCmd;
  delete  sexCmd;
  delete  bodypartCmd;
  delete  endCmd;
  delete  phantomDir;
  delete  bpDir;
}

void G4HumanPhantomMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
  if( command == cleanCmd )
    { 
      myUserPhantom->CleanPhantom();
    } 
  if( command == modelCmd )
    { 
      myUserPhantom->SetPhantomModel(newValue); 
    }       
  if( command == sexCmd )
    { 
      myUserPhantom->SetPhantomSex(newValue); 
    }
  if( command == bodypartCmd )
    {
      AddBodyPart(newValue);
    }
  if( command == endCmd )
    { 
      G4cout << 
	" ****************>>>> NEW PHANTOM CONSTRUCTION <<<<***************** " 
	     << G4endl;

      //myUserPhantom->UpdatePhantom();

      G4RunManager::GetRunManager()->Initialize();
    }
}

void  G4HumanPhantomMessenger::AddBodyPart(G4String newBodyPartSensitivity)
{

  const char* str[] = {newBodyPartSensitivity,NULL};

  std::string bodypart = strtok(str[0]," ");

  std::string sensitivity = strtok(NULL," ");

  if(sensitivity=="yes"){
    bps=true;
  }else{
    bps=false;
  }

  G4cout << " >>> Body Part = " << bodypart << "\n"
	 << " >>> Sensitivity = " << sensitivity << G4endl;


  myUserPhantom->SetBodyPartSensitivity(bodypart,bps);

}

