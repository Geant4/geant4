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
//
// Code developed by:
//  S.Guatelli, A. Le
//
//    *********************************
//    *                               *
//    *    BrachyDetectorMessenger.cc *
//    *                               *
//    *********************************
//
//
//
// 

#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"


BrachyDetectorMessenger::BrachyDetectorMessenger( BrachyDetectorConstruction* Det): detector(Det)
{ 
  detectorDir = new G4UIdirectory("/phantom/");
  detectorDir -> SetGuidance(" phantom control.");
      
  phantomMaterialCmd = new G4UIcmdWithAString("/phantom/selectMaterial",this);
  phantomMaterialCmd -> SetGuidance("Select Material of the phantom.");
  phantomMaterialCmd -> SetParameterName("choice",false);
  phantomMaterialCmd -> AvailableForStates(G4State_Idle);
  
  sourceCmd = new G4UIcmdWithAString("/source/switch",this);
  sourceCmd -> SetGuidance("Assign the selected geometry to G4RunManager."); 
  sourceCmd -> SetParameterName("choice",true);
  sourceCmd -> SetDefaultValue(" ");
  sourceCmd -> SetCandidates("TG186 Flexi Iodine Leipzig Oncura");
  sourceCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 
 }

BrachyDetectorMessenger::~BrachyDetectorMessenger()
{
  delete sourceCmd;
  delete phantomMaterialCmd; 
  delete detectorDir;
}

void BrachyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  // Change the material of the phantom
  if( command == phantomMaterialCmd )
   { detector -> SetPhantomMaterial(newValue);}

  // Switch the source in the phantom
  if( command == sourceCmd )
   {
    if(newValue=="Iodine" || newValue=="TG186"|| newValue=="Leipzig" || newValue== "Flexi" || newValue== "Oncura")
     { 
       detector -> SelectBrachytherapicSeed(newValue); 
       detector -> SwitchBrachytherapicSeed();
      }
   }
}

