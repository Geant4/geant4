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

#include "GammaKnifePhysicsListMessenger.hh"
#include "GammaKnifePhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

GammaKnifePhysicsListMessenger::GammaKnifePhysicsListMessenger(GammaKnifePhysicsList * physList)
:physicsList(physList)
{  
 listDir = new G4UIdirectory("/Physics/");
  // Building modular PhysicsList

 physicsListCmd = new G4UIcmdWithAString("/Physics/addPhysics",this);  
 physicsListCmd->SetGuidance("Add chunks of PhysicsList.");
 physicsListCmd->SetParameterName("physList",false);
 physicsListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
 gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setGCut",this);  
 gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setECut",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  protoCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setPCut",this);  
  protoCutCmd->SetGuidance("Set positron cut.");
  protoCutCmd->SetParameterName("Pcut",false);
  protoCutCmd->SetUnitCategory("Length");
  protoCutCmd->SetRange("Pcut>0.0");
  protoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/Physics/setCuts",this);  
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
  allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  


}

GammaKnifePhysicsListMessenger::~GammaKnifePhysicsListMessenger()
{
  delete physicsListCmd;
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
  delete allCutCmd;
  delete listDir;
}

void GammaKnifePhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if (command == physicsListCmd)
    { 
		physicsList->AddPhysicsList(newValue);
	}

 else if( command == gammaCutCmd )
   { physicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}
     
  else if( command == electCutCmd )
   { physicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));}
     
  else if( command == protoCutCmd )
   { physicsList->SetCutForPositron(protoCutCmd->GetNewDoubleValue(newValue));}

  else if( command == allCutCmd )
    {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      physicsList->SetCutForGamma(cut);
      physicsList->SetCutForElectron(cut);
      physicsList->SetCutForPositron(cut);
    } 

}






