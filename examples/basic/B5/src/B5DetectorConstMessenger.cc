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
/// \file B5DetectorConstMessenger.cc
/// \brief Implementation of the B5DetectorConstMessenger class

#include "B5DetectorConstMessenger.hh"
#include "B5DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstMessenger::B5DetectorConstMessenger(B5DetectorConstruction* det)
: G4UImessenger(), fDetector(det), fMydetDirectory(0), fArmCmd(0)
{
    fMydetDirectory = new G4UIdirectory("/mydet/");
    fMydetDirectory->SetGuidance("B5 detector setup control commands.");
    
    fArmCmd = new G4UIcmdWithADoubleAndUnit("/mydet/armAngle",this);
    fArmCmd->SetGuidance("Rotation angle of the second arm.");
    fArmCmd->SetParameterName("angle",true);
    fArmCmd->SetRange("angle>=0. && angle<180.");
    fArmCmd->SetDefaultValue(30.);
    fArmCmd->SetDefaultUnit("deg");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstMessenger::~B5DetectorConstMessenger()
{
    delete fArmCmd;
    delete fMydetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstMessenger::SetNewValue(G4UIcommand * command,
                                           G4String newValue)
{
    if( command==fArmCmd )
    { fDetector->SetArmAngle(fArmCmd->GetNewDoubleValue(newValue)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B5DetectorConstMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    if( command==fArmCmd )
    { cv = fArmCmd->ConvertToString(fDetector->GetArmAngle(),"deg"); }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
