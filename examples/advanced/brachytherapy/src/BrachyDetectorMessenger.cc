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
//
// Code developed by:
//  S.Guatelli
//
//    *********************************
//    *                               *
//    *    BrachyDetectorMessenger.cc *
//    *                               *
//    *********************************
//
//
// $Id: BrachyDetectorMessenger.cc,v 1.5 2002-12-06 16:32:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BrachyDetectorMessenger.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyRunAction.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BrachyDetectorMessenger::BrachyDetectorMessenger( BrachyDetectorConstruction* Det):
 Detector(Det)
{ 
  detDir = new G4UIdirectory("/detector/");
  detDir->SetGuidance(" detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/detector/setMaterial",this);
  AbsMaterCmd->SetGuidance("Select Material of the detector.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(G4State_Idle);
   selDetCmd = new G4UIcmdWithAString("/geom/select",this);
  
  mydetDir = new G4UIdirectory("/geom/");
  mydetDir->SetGuidance("Geometry setup commands.");
 
  selDetCmd->SetGuidance("Select the way detector geometry is built.");
  selDetCmd->SetGuidance(" Iodium:  Iodium Source ");
  selDetCmd->SetGuidance("  Iridium: Iridium  Source  "  );
  selDetCmd->SetGuidance("  Leipzig: Leipzig Source  "  );
   selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("Iridium");
  selDetCmd->SetCandidates("Iridium / Iodium / Leipzig");
  selDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  switchCmd = new G4UIcmdWithAString("/geom/switch",this);
  switchCmd->SetGuidance("Assign the selected geometry to G4RunManager.");
  switchCmd->SetGuidance("In case the choice is present to this command,");
  switchCmd->SetGuidance("\"/geom/select\" will be invoked and then switched.");
  switchCmd->SetParameterName("choice",true);
  switchCmd->SetDefaultValue(" ");
  switchCmd->SetCandidates("Iridium Iodium Leipzig ");
  switchCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BrachyDetectorMessenger::~BrachyDetectorMessenger()
{
 
  delete AbsMaterCmd; 
  delete   mydetDir;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BrachyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { Detector->SetAbsorberMaterial(newValue);}

 if( command == selDetCmd )
  {
    Detector->SelectDetector(newValue);
  }
  if( command == switchCmd )
  {
    if(newValue=="Iodium" || newValue=="Iridium"|| newValue=="Leipzig")
      { 
      Detector->SelectDetector(newValue); 
      Detector->SwitchDetector();
       }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
