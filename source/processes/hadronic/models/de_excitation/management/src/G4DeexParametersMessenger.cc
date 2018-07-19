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
// $Id: G4DeexParametersMessenger.cc 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4DeexParametersMessenger
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17-10-2017 
//
// Modifications:
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DeexParametersMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImanager.hh"
#include "G4DeexPrecoParameters.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DeexParametersMessenger::G4DeexParametersMessenger(G4DeexPrecoParameters* ptr) 
  : theParameters(ptr)
{
  fDirectory = new G4UIdirectory("/process/deex/");
  fDirectory->SetGuidance("Commands for nuclear de-excitation module.");

  readCmd = new G4UIcmdWithABool("/process/deex/readICdata",this);
  readCmd->SetGuidance("Enable/disable download IC data per atomic shell.");
  readCmd->SetParameterName("readIC",true);
  readCmd->SetDefaultValue(false);
  readCmd->AvailableForStates(G4State_PreInit);

  icCmd = new G4UIcmdWithABool("/process/deex/setIC",this);
  icCmd->SetGuidance("Enable/disable simulation of e- internal conversion.");
  icCmd->SetParameterName("IC",true);
  icCmd->SetDefaultValue(true);
  icCmd->AvailableForStates(G4State_PreInit);

  corgCmd = new G4UIcmdWithABool("/process/deex/correlatedGamma",this);
  corgCmd->SetGuidance("Enable/disable simulation of correlated gamma emission.");
  corgCmd->SetParameterName("corrG",true);
  corgCmd->SetDefaultValue(false);
  corgCmd->AvailableForStates(G4State_PreInit);

  maxjCmd = new G4UIcmdWithAnInteger("/process/deex/maxTwoJ",this);
  maxjCmd->SetGuidance("Set max value for 2J for simulation of correlated gamma emission.");
  maxjCmd->SetParameterName("max2J",true);
  maxjCmd->SetDefaultValue(10);
  maxjCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DeexParametersMessenger::~G4DeexParametersMessenger()
{
  delete fDirectory;

  delete readCmd;
  delete icCmd;
  delete corgCmd;
  delete maxjCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DeexParametersMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{
  if (command == readCmd) {
    theParameters->SetStoreICLevelData(readCmd->GetNewBoolValue(newValue));
  } else if (command == icCmd) {
    theParameters->SetInternalConversionFlag(icCmd->GetNewBoolValue(newValue));
  } else if (command == corgCmd) {
    theParameters->SetCorrelatedGamma(corgCmd->GetNewBoolValue(newValue));
  } else if (command == maxjCmd) { 
    theParameters->SetTwoJMAX(maxjCmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
