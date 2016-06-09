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
// $Id: BrachyDetectorMessenger.cc,v 1.8 2003/05/26 09:20:14 guatelli Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

#include "BrachyDetectorMessenger.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyRunAction.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

BrachyDetectorMessenger::BrachyDetectorMessenger( BrachyDetectorConstruction* Det): detector(Det)
{ 
  detectorDir = new G4UIdirectory("/phantom/");
  detectorDir->SetGuidance(" phantom control.");
      
  phantomMaterialCmd = new G4UIcmdWithAString("/phantom/selectMaterial",this);
  phantomMaterialCmd->SetGuidance("Select Material of the detector.");
  phantomMaterialCmd->SetParameterName("choice",false);
  phantomMaterialCmd->AvailableForStates(G4State_Idle);
  
  sourceCmd = new G4UIcmdWithAString("/source/switch",this);
  sourceCmd->SetGuidance("Assign the selected geometry to G4RunManager."); 
  sourceCmd->SetParameterName("choice",true);
  sourceCmd->SetDefaultValue(" ");
  sourceCmd->SetCandidates("Iridium Iodium Leipzig ");
  sourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 }

BrachyDetectorMessenger::~BrachyDetectorMessenger()
{
  delete sourceCmd;
  delete phantomMaterialCmd; 
  delete detectorDir;
}

void BrachyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == phantomMaterialCmd )
   { detector->SetPhantomMaterial(newValue);}

  if( command == sourceCmd )
   {
    if(newValue=="Iodium" || newValue=="Iridium"|| newValue=="Leipzig")
     { 
       detector->SelectBrachytherapicSeed(newValue); 
       detector->SwitchBrachytherapicSeed();
      }
   }
}

