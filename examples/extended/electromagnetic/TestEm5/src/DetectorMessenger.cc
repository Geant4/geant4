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
/// \file electromagnetic/TestEm5/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),
 fAbsMaterCmd(0),
 fAbsThickCmd(0),
 fAbsSizYZCmd(0),
 fAbsXposCmd(0),
 fWorldMaterCmd(0),
 fWorldXCmd(0),
 fWorldYZCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("UI commands specific to this example.");
  
  fDetDir = new G4UIdirectory("/testem/det/", false);
  fDetDir->SetGuidance("detector construction commands");
      
  fAbsMaterCmd = new G4UIcmdWithAString("/testem/det/setAbsMat",this);
  fAbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  fAbsMaterCmd->SetParameterName("choice",false);
  fAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsMaterCmd->SetToBeBroadcasted(false);
  
  fWorldMaterCmd = new G4UIcmdWithAString("/testem/det/setWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select Material of the World.");
  fWorldMaterCmd->SetParameterName("wchoice",false);
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldMaterCmd->SetToBeBroadcasted(false);
    
  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsThick",this);
  fAbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  fAbsThickCmd->SetParameterName("SizeZ",false);  
  fAbsThickCmd->SetRange("SizeZ>0.");
  fAbsThickCmd->SetUnitCategory("Length");  
  fAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsThickCmd->SetToBeBroadcasted(false);
  
  fAbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsYZ",this);
  fAbsSizYZCmd->SetGuidance("Set sizeYZ of the Absorber");
  fAbsSizYZCmd->SetParameterName("SizeYZ",false);
  fAbsSizYZCmd->SetRange("SizeYZ>0.");
  fAbsSizYZCmd->SetUnitCategory("Length");
  fAbsSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsSizYZCmd->SetToBeBroadcasted(false);
  
  fAbsXposCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsXpos",this);
  fAbsXposCmd->SetGuidance("Set X pos. of the Absorber");
  fAbsXposCmd->SetParameterName("Xpos",false);
  fAbsXposCmd->SetUnitCategory("Length");
  fAbsXposCmd->AvailableForStates(G4State_PreInit);
  fAbsXposCmd->SetToBeBroadcasted(false);
  
  fWorldXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldX",this);
  fWorldXCmd->SetGuidance("Set X size of the World");
  fWorldXCmd->SetParameterName("WSizeX",false);
  fWorldXCmd->SetRange("WSizeX>0.");
  fWorldXCmd->SetUnitCategory("Length");
  fWorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldXCmd->SetToBeBroadcasted(false);
  
  fWorldYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldYZ",this);
  fWorldYZCmd->SetGuidance("Set sizeYZ of the World");
  fWorldYZCmd->SetParameterName("WSizeYZ",false);
  fWorldYZCmd->SetRange("WSizeYZ>0.");
  fWorldYZCmd->SetUnitCategory("Length");
  fWorldYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldYZCmd->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fAbsMaterCmd; 
  delete fAbsThickCmd; 
  delete fAbsSizYZCmd;  
  delete fAbsXposCmd; 
  delete fWorldMaterCmd;
  delete fWorldXCmd;
  delete fWorldYZCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == fAbsMaterCmd )
   {fDetector->SetAbsorberMaterial(newValue);}
   
  if ( command == fWorldMaterCmd )
   {fDetector->SetWorldMaterial(newValue);}
   
  if ( command == fAbsThickCmd )
  {fDetector->SetAbsorberThickness(fAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fAbsSizYZCmd )
   {fDetector->SetAbsorberSizeYZ(fAbsSizYZCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fAbsXposCmd )
   {fDetector->SetAbsorberXpos(fAbsXposCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fWorldXCmd )
   {fDetector->SetWorldSizeX(fWorldXCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fWorldYZCmd )
   {fDetector->SetWorldSizeYZ(fWorldYZCmd->GetNewDoubleValue(newValue));}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
