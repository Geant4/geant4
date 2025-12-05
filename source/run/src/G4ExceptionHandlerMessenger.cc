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
// G4ExceptionHandlerMessenger implementation
//
// Original author: M.Asai, 1997
// --------------------------------------------------------------------

#include "G4ExceptionHandlerMessenger.hh"

#include "G4ExceptionHandler.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UIparameter.hh"

G4ExceptionHandlerMessenger::G4ExceptionHandlerMessenger(G4ExceptionHandler* expH) 
: expHandler(expH)
{
  EHDirectory = new G4UIdirectory("/control/exception/");
  EHDirectory->SetGuidance("Managing number of G4Exception warning messages");

  nwCmd = new G4UIcommand("/control/exception/maxExceptionWarning",this);
  nwCmd->SetGuidance("Set maximum number of G4Exception warning messages to be displayed");
  nwCmd->SetGuidance("for the specified error code.");
  nwCmd->SetGuidance("Number of warnings is counted for each thread individually.");
  nwCmd->SetGuidance("Once the number reaches to the maximum, the warning message of specified");
  nwCmd->SetGuidance("error code won't be printed out any more.");
  nwCmd->SetGuidance("Number is counted through the program execution if more than one run");
  nwCmd->SetGuidance("are executed. Count can be reset with the same command.");
  nwCmd->SetGuidance("If number is set to zero, warning message of the specified error code");
  nwCmd->SetGuidance("won't be displayed at all.");
  nwCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  auto p1 = new G4UIparameter("errorCode",'s',false);
  nwCmd->SetParameter(p1);
  auto p2 = new G4UIparameter("maxNumber",'i',false);
  p2->SetParameterRange("maxNumber >= 0");
  nwCmd->SetParameter(p2);

  tnwCmd = new G4UIcmdWithAnInteger("/control/exception/maxTotalExceptionWarning",this);
  tnwCmd->SetGuidance("Set total maximum number of G4Exception warning messages to be displayed.");
  tnwCmd->SetGuidance("Number of warnings is counted for each thread individually.");
  tnwCmd->SetGuidance("Once the number reaches to the maximum, warning messages won't be"); 
  tnwCmd->SetGuidance("printed out any more.");
  tnwCmd->SetGuidance("Number is counted through the program execution if more than one run");
  tnwCmd->SetGuidance("are executed. Count can be reset with the same command.");
  tnwCmd->SetGuidance("If number is set to zero, warning message won't be displayed at all.");
  tnwCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  tnwCmd->SetParameterName("maxNumber",false);
  tnwCmd->SetRange("maxNumber >= 0");

}

// --------------------------------------------------------------------
G4ExceptionHandlerMessenger::~G4ExceptionHandlerMessenger()
{
  delete nwCmd;
  delete tnwCmd;
  delete EHDirectory;
}

// --------------------------------------------------------------------
void G4ExceptionHandlerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == nwCmd) {
    G4Tokenizer next(newValue);
    G4String ec = next();
    G4int mc = StoI(next());
    expHandler->SetMaxWarning(ec,mc);
  }
  else if (command == tnwCmd) {
    expHandler->SetMaxTotalWarning(tnwCmd->GetNewIntValue(newValue));
  }
}

// --------------------------------------------------------------------
G4String G4ExceptionHandlerMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  G4String cv;
  return cv;
}

