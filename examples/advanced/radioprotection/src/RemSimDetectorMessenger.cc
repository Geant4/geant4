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
//    *    RemSimDetectorMessenger.cc *
//    *                               *
//    *********************************
//
//
// $Id: RemSimDetectorMessenger.cc,v 1.1 2004-01-30 12:25:44 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "RemSimDetectorMessenger.hh"
#include "RemSimRunAction.hh"
#include "RemSimDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

RemSimDetectorMessenger::RemSimDetectorMessenger( RemSimDetectorConstruction* Det): detector(Det)
{ 
  vehicleDir = new G4UIdirectory("/vehicle/");
  vehicleDir->SetGuidance("vehicle control.");
        
  vehicleCmd = new G4UIcmdWithAString("/vehicle/switch",this);
  vehicleCmd->SetGuidance("Assign the selected vehicle to G4RunManager."); 
  vehicleCmd->SetParameterName("choice",true);
  vehicleCmd->SetCandidates("Vehicle1 Vehicle2 ");
  vehicleCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 }

RemSimDetectorMessenger::~RemSimDetectorMessenger()
{
  delete vehicleCmd;
  delete vehicleDir;
}

void RemSimDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   if( command == vehicleCmd )
   {
     if(newValue=="Vehicle1" || newValue=="Vehicle2") 
       detector->SwitchVehicle(newValue);
     else G4cout<< "This vehicle concept is not available!"<<G4cout; 
  }
}

