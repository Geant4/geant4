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
// $Id: WLSEventActionMessenger.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/WLSEventActionMessenger.cc
/// \brief Implementation of the WLSEventActionMessenger class
//
//
#include "globals.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "WLSEventAction.hh"
#include "WLSEventActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventActionMessenger::WLSEventActionMessenger(WLSEventAction* eventaction)
  : fEventAction(eventaction)
{
  fSetVerboseCmd = new G4UIcmdWithAnInteger("/event/setverbose",this);
  fSetVerboseCmd->SetGuidance("Set verbose level ." );
  fSetVerboseCmd->SetParameterName("level",true);
  fSetVerboseCmd->SetDefaultValue(0);

  fDrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  fDrawCmd->SetGuidance("Draw the tracks in the event");
  fDrawCmd->SetGuidance("  Choice : none,charged, all");
  fDrawCmd->SetParameterName("choice",true);
  fDrawCmd->SetDefaultValue("all");
  fDrawCmd->SetCandidates("none charged all");
  fDrawCmd->AvailableForStates(G4State_Idle);

  fPrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  fPrintCmd->SetGuidance("Print events modulo n");
  fPrintCmd->SetParameterName("EventNb",false);
  fPrintCmd->SetRange("EventNb>0");
  fPrintCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventActionMessenger::~WLSEventActionMessenger()
{
  delete fSetVerboseCmd;
  delete fDrawCmd;
  delete fPrintCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if (command == fSetVerboseCmd)
    fEventAction->SetEventVerbose(fSetVerboseCmd->GetNewIntValue(newValue));
 
  if (command == fDrawCmd)
    fEventAction->SetDrawFlag(newValue);

  if (command == fPrintCmd)
    fEventAction->SetPrintModulo(fPrintCmd->GetNewIntValue(newValue));
}
