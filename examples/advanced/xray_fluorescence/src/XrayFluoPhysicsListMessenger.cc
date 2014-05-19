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
// History:
// 18 Jan 2011 Alf Adapted from TestEm18
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "XrayFluoPhysicsListMessenger.hh"

#include "XrayFluoPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayFluoPhysicsListMessenger::XrayFluoPhysicsListMessenger(XrayFluoPhysicsList* pPhys)
:pPhysicsList(pPhys)
{
  physDir = new G4UIdirectory("/phys/");
  physDir->SetGuidance("physics list commands");
  
  pListCmd = new G4UIcmdWithAString("/phys/addPhysics",this);  
  pListCmd->SetGuidance("Add modules physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit);
  
  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/phys/setGCut",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/phys/setECut",this);  
  electCutCmd->SetGuidance("Set electron and positron cuts.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  protoCutCmd = new G4UIcmdWithADoubleAndUnit("/phys/setPCut",this);  
  protoCutCmd->SetGuidance("Set proton cut.");
  protoCutCmd->SetParameterName("Pcut",false);
  protoCutCmd->SetUnitCategory("Length");
  protoCutCmd->SetRange("Pcut>0.0");
  protoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  allCutCmd = new G4UIcmdWithADoubleAndUnit("/phys/setCuts",this);  
  allCutCmd->SetGuidance("Set cut for all.");
  allCutCmd->SetParameterName("cut",false);
  allCutCmd->SetUnitCategory("Length");
  allCutCmd->SetRange("cut>0.0");
  allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  /*fluoCmd = new G4UIcmdWithABool("/phys/fluo",this);  
  fluoCmd->SetGuidance("Set fluorescence on/off.");
  fluoCmd->SetParameterName("fluo",false);
  fluoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  pixeCmd = new G4UIcmdWithABool("/phys/pixe",this);  
  pixeCmd->SetGuidance("Set PIXE on/off.");
  pixeCmd->SetParameterName("pixe",false);
  pixeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayFluoPhysicsListMessenger::~XrayFluoPhysicsListMessenger()
{
  delete pListCmd;
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
  delete allCutCmd;
  delete physDir;
  //  delete fluoCmd;
  //  delete pixeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayFluoPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
         
  if( command == gammaCutCmd )
   { pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == electCutCmd )
   { 
     G4double cut = electCutCmd->GetNewDoubleValue(newValue);
     pPhysicsList->SetCutForElectron(cut);
     pPhysicsList->SetCutForPositron(cut);
   }
     
  if( command == protoCutCmd )
   { pPhysicsList->SetCutForProton(protoCutCmd->GetNewDoubleValue(newValue));}

  if( command == allCutCmd )
    {
      G4double cut = allCutCmd->GetNewDoubleValue(newValue);
      pPhysicsList->SetCutForGamma(cut);
      pPhysicsList->SetCutForElectron(cut);
      pPhysicsList->SetCutForPositron(cut);
      pPhysicsList->SetCutForProton(cut);
    } 

  //Notify the run manager that the physics has been modified
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

//  if( command == fluoCmd )
//   { pPhysicsList->SetFluorescence(fluoCmd->GetNewBoolValue(newValue));}
//
//  if( command == pixeCmd )
//   { pPhysicsList->SetPIXE(fluoCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
