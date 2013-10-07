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
/// \file B5MagneticFieldMessenger.cc
/// \brief Implementation of the B5MagneticFieldMessenger class

#include "B5MagneticFieldMessenger.hh"
#include "B5MagneticField.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5MagneticFieldMessenger::B5MagneticFieldMessenger(B5MagneticField* field)
: G4UImessenger(), fField(field), fFieldCmd(0)
{
    fFieldCmd = new G4UIcmdWithADoubleAndUnit("/mydet/fieldValue",this);
    fFieldCmd->SetGuidance("Field strength");
    fFieldCmd->SetParameterName("field",true);
    fFieldCmd->SetDefaultValue(1.);
    fFieldCmd->SetDefaultUnit("tesla");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5MagneticFieldMessenger::~B5MagneticFieldMessenger()
{
    delete fFieldCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5MagneticFieldMessenger::SetNewValue(G4UIcommand * command,
                                           G4String newValue)
{
    if( command==fFieldCmd )
    { fField->SetField(fFieldCmd->GetNewDoubleValue(newValue)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String B5MagneticFieldMessenger::GetCurrentValue(G4UIcommand * command)
{
    G4String cv;
    if( command==fFieldCmd )
    { cv = fFieldCmd->ConvertToString(fField->GetField(),"tesla"); }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
