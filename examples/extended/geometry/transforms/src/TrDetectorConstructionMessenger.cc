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
// $Id$
//
/// \file TrDetectorConstructionMessenger.cc
/// \brief Implementation of the TrDetectorConstructionMessenger class

#include "TrDetectorConstructionMessenger.hh"
#include "TrDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrDetectorConstructionMessenger::TrDetectorConstructionMessenger(
                           TrDetectorConstruction* detectorConstruction)
 : G4UImessenger(),
   fDetectorConstruction(detectorConstruction),
   fDirectory(0),
   fSetMethodCmd(0)
{ 
  fDirectory = new G4UIdirectory("/transform/");
  fDirectory->SetGuidance("Transform example detector control");
       
  fSetMethodCmd 
    = new G4UIcmdWithAString("/transform/setMethod",this);
  fSetMethodCmd->SetGuidance("Select method for definition of transformations.");
  fSetMethodCmd->SetParameterName("Method", false);
  fSetMethodCmd->SetCandidates(
    "WithDirectMatrix WithInverseMatrix WithAxialRotations WithEulerAngles WithReflections");
  fSetMethodCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrDetectorConstructionMessenger::~TrDetectorConstructionMessenger()
{
  delete fDirectory;
  delete fSetMethodCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrDetectorConstructionMessenger::SetNewValue(G4UIcommand* command, 
                                                     G4String newValue)
{ 
  if( command == fSetMethodCmd ) { 
    TrDetectorConstruction::EMethod method;
    if ( newValue == "WithDirectMatrix" )         
      method = TrDetectorConstruction::kWithDirectMatrix;
    else if ( newValue == "WithInverseMatrix" )         
      method = TrDetectorConstruction::kWithInverseMatrix;
    else if ( newValue == "WithAxialRotations" )     
      method = TrDetectorConstruction::kWithAxialRotations;
    else if ( newValue == "WithEulerAngles" ) 
      method = TrDetectorConstruction::kWithEulerAngles;
    else if ( newValue == "WithReflections" )    
      method = TrDetectorConstruction::kWithReflections;
    else return;  
    fDetectorConstruction->SetMethod(method);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
