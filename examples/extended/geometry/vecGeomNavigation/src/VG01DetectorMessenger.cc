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
/// \file VG01DetectorMessenger.cc
/// \brief Implementation of the VG01DetectorMessenger class

//  Authors: J. Apostolakis & S. Wenzel (CERN)  2018-2021
//
//  Started from FullCMS code by Mihaly Novak (CERN) 2017  

#include "VG01DetectorMessenger.hh"

#include "VG01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"


VG01DetectorMessenger::
VG01DetectorMessenger( VG01DetectorConstruction* myDet )
  : G4UImessenger(), fDetector( myDet ) 
{

  fDetectorDir = new G4UIdirectory( "/mydet/" );
  fDetectorDir->SetGuidance( "Detector control." );

  fFieldCommand = new G4UIcmdWithADoubleAndUnit( "/mydet/setField", this );
  fFieldCommand->SetGuidance( "Define uniform magnetic field along Z." );
  fFieldCommand->SetGuidance( " -> in unit of  [Tesla]" );
  fFieldCommand->SetParameterName( "By", false );
  fFieldCommand->SetDefaultValue( 0.0 );
  fFieldCommand->SetUnitCategory( "Magnetic flux density" );
  fFieldCommand->AvailableForStates( G4State_PreInit, G4State_Idle );

  fGDMLCommand = new G4UIcmdWithAString( "/mydet/setGdmlFile", this );
  fGDMLCommand->SetGuidance( "Set the GDML file." );
  fGDMLCommand->SetDefaultValue( "TestNTST.gdml" );
  fGDMLCommand->AvailableForStates( G4State_PreInit, G4State_Idle );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VG01DetectorMessenger::~VG01DetectorMessenger() {
  delete fFieldCommand;
  delete fDetectorDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) 
{
  if ( command == fFieldCommand ) {
    fDetector->SetMagFieldValue(fFieldCommand->GetNewDoubleValue(newValue));
  }
  else
  {
    if ( command == fGDMLCommand ) {
      fDetector->SetGDMLFileName( newValue );
    } else {
      G4cerr << "VG01DetectorMessenger: ERROR> Unknown command " << G4endl;
    }
  } 
}
