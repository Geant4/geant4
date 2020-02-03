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
//  B. Tomé
//
//    *********************************
//    *                               *
//    *    UltraDetectorMessenger.cc *
//    *                               *
//    *********************************
//
//
//
// 

#include "UltraDetectorMessenger.hh"
#include "UltraDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"


UltraDetectorMessenger::UltraDetectorMessenger( UltraDetectorConstruction* Det): fDetector(Det)
{ 
  fDetectorDir = new G4UIdirectory("/detector/");
  fDetectorDir -> SetGuidance("Configuration of Ultra setup.");
  
  fMirrorCmd = new G4UIcmdWithAString("/detector/reflection",this);
  fMirrorCmd -> SetGuidance("Define the type of reflecting surface. Default is none."); 
  fMirrorCmd -> SetParameterName("choice",true);
  fMirrorCmd -> SetDefaultValue("none");
  fMirrorCmd -> SetCandidates("none specular ground");
  fMirrorCmd -> AvailableForStates(G4State_Idle); 
 }

UltraDetectorMessenger::~UltraDetectorMessenger()
{
  delete fMirrorCmd;
  delete fDetectorDir;
}

void UltraDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  // Select the type of reflection at ground
  if( command == fMirrorCmd ) {
    fDetector -> SetReflectionType(newValue); 
  }
}

