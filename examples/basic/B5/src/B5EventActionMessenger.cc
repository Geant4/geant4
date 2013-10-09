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
/// \file B5EventActionMessenger.cc
/// \brief Implementation of the B5EventActionMessenger class

#include "B5EventActionMessenger.hh"
#include "B5EventAction.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventActionMessenger::B5EventActionMessenger(B5EventAction* mpga)
: G4UImessenger(), fTarget(mpga), fVerboseCmd(0)
{
    fVerboseCmd = new G4UIcmdWithAnInteger("/mydet/verbose",this);
    fVerboseCmd->SetGuidance("Event verbose level.");
    fVerboseCmd->SetGuidance(
                " Event summary will be displayed for every 'level' events.");
    fVerboseCmd->SetParameterName("level",true);
    fVerboseCmd->SetRange("level>=0");
    fVerboseCmd->SetDefaultValue(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventActionMessenger::~B5EventActionMessenger()
{
    delete fVerboseCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventActionMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue)
{
    if ( command==fVerboseCmd )
    { fTarget ->SetVerbose(fVerboseCmd->GetNewIntValue(newValue)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B5EventActionMessenger::GetCurrentValue(G4UIcommand* command)
{
    G4String cv;
    if ( command==fVerboseCmd )
    { cv = fVerboseCmd->ConvertToString(fTarget ->GetVerbose()); }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
