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
// Previous authors: G. Guerrieri, S. Guatelli, and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
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

  modelCmd = new G4UIcmdWithAString("/phantom/setPhantomModel",this);
  modelCmd->SetGuidance("Set sex of Phantom: MIRD, ORNLFemale, ORNLMale, MIX, MIRDHead, ORNLHead.");
  modelCmd->SetParameterName("phantomModel",true);
  modelCmd->SetDefaultValue("MIRD");
  modelCmd->SetCandidates("MIRD ORNLFemale ORNLMale MIX MIRDHead ORNLHead");
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
  delete  modelCmd;
  delete  sexCmd;
  delete  bodypartCmd;
  delete  endCmd;
  delete  phantomDir;
  delete  bpDir;
}

void G4HumanPhantomMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 

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
    }
}

void  G4HumanPhantomMessenger::AddBodyPart(G4String newBodyPartSensitivity)
{

  char* str = new char[newBodyPartSensitivity.length()+1];

  strcpy(str, newBodyPartSensitivity.c_str()); 
  
  std::string bpart = strtok(str," ");

  std::string sensitivity = strtok(NULL," ");

  if(sensitivity=="yes"){
    bps=true;
  }else{
    bps=false;
  }

  G4cout << " >>> Body Part = " << bpart << "\n"
	 << " >>> Sensitivity = " << sensitivity << G4endl;

  myUserPhantom->SetBodyPartSensitivity(bpart,bps);
}

