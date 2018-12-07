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
/// \file persistency/gdml/G02/src/G02DetectorMessenger.cc
/// \brief Implementation of the G02DetectorMessenger class
//
//
//
// Class G02DetectorMessenger implementation
//
// ----------------------------------------------------------------------------

#include "globals.hh"

#include "G02DetectorMessenger.hh"
#include "G02DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02DetectorMessenger::G02DetectorMessenger( G02DetectorConstruction* myDet )
  : G4UImessenger(),
    fTheDetector( myDet ),
    fTheDetectorDir(0),
    fTheReadCommand(0),
    fTheWriteCommand(0),
    fTheStepCommand(0)
{ 
  fTheDetectorDir = new G4UIdirectory( "/mydet/" );
  fTheDetectorDir->SetGuidance("Detector control.");

  fTheReadCommand = new G4UIcmdWithAString("/mydet/readFile", this);
  fTheReadCommand ->SetGuidance("READ GDML file with given name");
  fTheReadCommand ->SetParameterName("FileRead", false);
  fTheReadCommand ->SetDefaultValue("test.gdml");
  fTheReadCommand ->AvailableForStates(G4State_PreInit);
  
  fTheWriteCommand = new G4UIcmdWithAString("/mydet/writeFile", this);
  fTheWriteCommand ->SetGuidance("WRITE geometry to GDML file with given name");
  fTheWriteCommand ->SetParameterName("FileWrite", false);
  fTheWriteCommand ->SetDefaultValue("wtest.gdml");
  fTheWriteCommand ->AvailableForStates(G4State_PreInit);

  fTheStepCommand = new G4UIcmdWithAString("/mydet/StepFile", this);
  fTheStepCommand ->SetGuidance("Read STEP Tools files (name without extension)");
  fTheStepCommand ->SetParameterName("STEPFile", false);
  fTheStepCommand ->SetDefaultValue("mbb");
  fTheStepCommand ->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02DetectorMessenger::~G02DetectorMessenger()
{
  delete fTheReadCommand;
  delete fTheWriteCommand;
  delete fTheStepCommand;
  delete fTheDetectorDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G02DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if ( command == fTheReadCommand )
  { 
    fTheDetector->SetReadFile(newValue );
  }
  if ( command == fTheWriteCommand )
  { 
    fTheDetector->SetWriteFile(newValue );
  }
  if ( command == fTheStepCommand )
  { 
    fTheDetector->SetStepFile(newValue );
  }
}
