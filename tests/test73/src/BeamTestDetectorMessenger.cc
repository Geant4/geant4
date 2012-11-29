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
//  BeamTestDetectorMessenger.cc
//  MSCTest
//
//  Created by Andrea Dotti on 3/7/12.
//

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "BeamTestDetectorMessenger.hh"
#include "BeamTestDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "BeamTestEventAction.hh"
#include "G4RunManager.hh"

BeamTestDetectorMessenger::BeamTestDetectorMessenger( BeamTestDetectorConstruction* detector) : 
    theDetector(detector),
    theDetectorDir(new G4UIdirectory("/mydet/")),
    theChamberWidthCmd(new G4UIcmdWithADoubleAndUnit("/mydet/chamberWidth",this)),
    theNumberOfChambersCmd(new G4UIcmdWithAnInteger("/mydet/numberOfChambers",this)),
    theChamberSpacingCmd(new G4UIcmdWithADoubleAndUnit("/mydet/chamberSpacing", this)),
    theUpdateCommand(new G4UIcmdWithoutParameter("/mydet/update",this))
{
    theDetectorDir->SetGuidance("Commands to control detector setup");
    
    theChamberWidthCmd->SetGuidance("Set chamber width (in um)");
    theChamberWidthCmd->SetDefaultValue(1*CLHEP::um);
    theChamberWidthCmd->SetDefaultUnit("um");
    theChamberWidthCmd->SetParameterName("CW", false);
    theChamberWidthCmd->SetUnitCategory("Length");
    theChamberWidthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    theChamberSpacingCmd->SetGuidance("Set chambers spacing (in um)");
    theChamberSpacingCmd->SetDefaultValue(1*CLHEP::um);
    theChamberSpacingCmd->SetDefaultUnit("um");
    theChamberSpacingCmd->SetParameterName("CS", false);
    theChamberSpacingCmd->SetUnitCategory("Length");
    theChamberSpacingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    theNumberOfChambersCmd->SetGuidance("Set number of chambers");
    theNumberOfChambersCmd->SetDefaultValue(1);
    theNumberOfChambersCmd->SetParameterName("NC", false);
    theNumberOfChambersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    theUpdateCommand->SetGuidance("Update geometry");
    theUpdateCommand->SetGuidance("This command MUST be applied before \"beamOn\" ");
    theUpdateCommand->SetGuidance("if you changed geometrical value(s).");
    theUpdateCommand->AvailableForStates(G4State_Idle);    
}
                   

BeamTestDetectorMessenger::~BeamTestDetectorMessenger() 
{ 
    delete theDetectorDir;
    delete theChamberSpacingCmd;
    delete theChamberWidthCmd;
    delete theNumberOfChambersCmd;
    delete theUpdateCommand;
}

void BeamTestDetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
    if ( command == theChamberSpacingCmd )
        theDetector->SetChamberSpacing(theChamberSpacingCmd->GetNewDoubleValue(newValue));
    
    if ( command == theChamberWidthCmd )
        theDetector->SetChamberWidth(theChamberWidthCmd->GetNewDoubleValue(newValue));
    
    if ( command == theNumberOfChambersCmd ) {
        theDetector->SetNumberChambers(theNumberOfChambersCmd->GetNewIntValue(newValue));
        //Need to update the analysis part with number of chambers..
        //This can be skipped if hits contain information on number of chamber itself...
        G4UserEventAction* theEventAction = const_cast<G4UserEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
        static_cast<BeamTestEventAction*>(theEventAction)->SetNumberOfChambers(theNumberOfChambersCmd->GetNewIntValue(newValue));
    }
    if ( command == theUpdateCommand )
        theDetector->UpdateGeometry();
}
