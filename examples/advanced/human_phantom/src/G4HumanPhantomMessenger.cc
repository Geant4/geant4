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
  :fUserPhantom(myUsrPhtm),fBps(false)
{ 
  fPhantomDir = new G4UIdirectory("/phantom/");
  fPhantomDir->SetGuidance("Set Your Phantom.");
  
  fDir = new G4UIdirectory("/bodypart/");
  fDir->SetGuidance("Add Body Part to Phantom");

  fModelCmd = new G4UIcmdWithAString("/phantom/setPhantomModel",this);
  fModelCmd->SetGuidance("Set sex of Phantom: MIRD, ORNLFemale, ORNLMale, MIX, MIRDHead, ORNLHead.");
  fModelCmd->SetParameterName("phantomModel",true);
  fModelCmd->SetDefaultValue("MIRD");
  fModelCmd->SetCandidates("MIRD ORNLFemale ORNLMale MIRDHead ORNLHead");
  fModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fSexCmd = new G4UIcmdWithAString("/phantom/setPhantomSex",this);
  fSexCmd->SetGuidance("Set sex of Phantom: Male or Female.");
  fSexCmd->SetParameterName("phantomSex",true);
  fSexCmd->SetDefaultValue("Female");
  fSexCmd->SetCandidates("Male Female");
  fSexCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  fBodypartCmd = new G4UIcmdWithAString("/bodypart/addBodyPart",this);
  fBodypartCmd->SetGuidance("Add a Body Part to Phantom");
  fBodypartCmd->SetParameterName("bpName",true);
  fBodypartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fEndCmd = new G4UIcmdWithoutParameter("/phantom/buildNewPhantom",this);
  fEndCmd->SetGuidance("Build your Phantom.");
  fEndCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

G4HumanPhantomMessenger::~G4HumanPhantomMessenger()
{
  delete  fModelCmd;
  delete  fSexCmd;
  delete  fBodypartCmd;
  delete  fEndCmd;
  delete  fPhantomDir;
  delete  fDir;
}

void G4HumanPhantomMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 

  if( command == fModelCmd )
    { 
      fUserPhantom->SetPhantomModel(newValue); 
    }       
  if( command == fSexCmd )
    { 
      fUserPhantom->SetPhantomSex(newValue); 
    }
  if( command == fBodypartCmd )
    {
      AddBodyPart(newValue);
    }
  if( command == fEndCmd )
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

  std::string sensitivity = strtok(nullptr," ");

  if(sensitivity=="yes"){
    fBps=true;
  }else{
    fBps=false;
  }

  G4cout << " >>> Body Part = " << bpart << "\n"
	 << " >>> Sensitivity = " << sensitivity << G4endl;

  fUserPhantom->SetBodyPartSensitivity(bpart,fBps);
}

