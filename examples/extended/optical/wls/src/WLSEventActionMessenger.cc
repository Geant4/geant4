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
/// \file optical/wls/src/WLSEventActionMessenger.cc
/// \brief Implementation of the WLSEventActionMessenger class
//
//
#include "globals.hh"

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventActionMessenger::~WLSEventActionMessenger()
{
  delete fSetVerboseCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if (command == fSetVerboseCmd)
    fEventAction->SetEventVerbose(fSetVerboseCmd->GetNewIntValue(newValue));
}
