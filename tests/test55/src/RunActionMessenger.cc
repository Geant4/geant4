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

#include "RunActionMessenger.hh"
#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include <sstream>


RunActionMessenger::RunActionMessenger(RunAction* run) :
   runAction(run) {

  runDirectory = new G4UIdirectory("/testem/run/");
  runDirectory -> SetGuidance("Run control");
 
  rangeTestCmd = new G4UIcommand("/testem/run/testRange",this);
  rangeTestCmd -> SetGuidance("Test simulated range against reference:");
  rangeTestCmd -> SetGuidance("  <ref.value> <abs.error>  <unit>");
  rangeTestCmd -> AvailableForStates(G4State_Idle);  

  G4UIparameter* refRange = new G4UIparameter("range", 'd', false);
  refRange -> SetGuidance("reference range");
  refRange -> SetParameterRange("range>0");
  rangeTestCmd -> SetParameter(refRange);

  G4UIparameter* refRangeErr = new G4UIparameter("rangeerror", 'd', false);
  refRangeErr -> SetGuidance("reference range abs. error");
  refRangeErr -> SetParameterRange("rangeerror>0");
  rangeTestCmd -> SetParameter(refRangeErr);

  G4UIparameter* refRangeUnit = new G4UIparameter("rangeunit", 's', false);
  refRangeUnit -> SetGuidance("reference range unit");
  rangeTestCmd -> SetParameter(refRangeUnit);

  elossTestCmd = new G4UIcommand("/testem/run/testMeanEnergyLoss",this);
  elossTestCmd -> SetGuidance("Test computed mean energy loss (of primary)");
  elossTestCmd -> SetGuidance("against reference:");
  elossTestCmd -> SetGuidance("  <ref.value> <abs.error>  <unit>");
  elossTestCmd -> AvailableForStates(G4State_Idle);  

  G4UIparameter* refEloss = new G4UIparameter("eloss", 'd', false);
  refEloss -> SetGuidance("reference mean energy loss");
  refEloss -> SetParameterRange("eloss>0");
  elossTestCmd -> SetParameter(refEloss);

  G4UIparameter* refElossErr = new G4UIparameter("elosserror", 'd', false);
  refElossErr -> SetGuidance("reference abs. error");
  refElossErr -> SetParameterRange("elosserror>0");
  elossTestCmd -> SetParameter(refElossErr);

  G4UIparameter* refElossUnit = new G4UIparameter("elossunit", 's', false);
  refElossUnit -> SetGuidance("reference unit");
  elossTestCmd -> SetParameter(refElossUnit);
}


RunActionMessenger::~RunActionMessenger() {

  delete rangeTestCmd;
  delete elossTestCmd;
  delete runDirectory;         
}


void RunActionMessenger::SetNewValue(G4UIcommand* command,
                                     G4String newValue) { 

  if( command == rangeTestCmd ) {

     std::istringstream isstream( newValue );
    
     G4double range;
     G4double error;
     G4String unit;

     isstream >> range >> error >> unit;
     G4double valUnit =  G4UIcommand::ValueOf(unit);
  
     runAction -> CreateRangeTest( range * valUnit, error/range );
  }
  if( command == elossTestCmd ) {

     std::istringstream isstream( newValue );
    
     G4double eloss;
     G4double error;
     G4String unit;

     isstream >> eloss >> error >> unit;
     G4double valUnit =  G4UIcommand::ValueOf(unit);
  
     runAction -> CreateEnergyLossTest( eloss * valUnit, error/eloss );
  }
}

