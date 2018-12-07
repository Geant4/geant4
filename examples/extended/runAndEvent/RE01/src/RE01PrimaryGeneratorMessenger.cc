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
/// \file runAndEvent/RE01/src/RE01PrimaryGeneratorMessenger.cc
/// \brief Implementation of the RE01PrimaryGeneratorMessenger class
//
//

#include "RE01PrimaryGeneratorMessenger.hh"
#include "RE01PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01PrimaryGeneratorMessenger
::RE01PrimaryGeneratorMessenger(RE01PrimaryGeneratorAction * mpga)
  :G4UImessenger(),fMyAction(mpga),
   fMydetDirectory(0),fGenCmd(0)
{
  fMydetDirectory = new G4UIdirectory("/mydet/");
  fMydetDirectory->SetGuidance("RE01 detector control commands.");

  fGenCmd = new G4UIcmdWithAString("/mydet/generator",this);
  fGenCmd->SetGuidance("Select primary generator.");
  fGenCmd->SetGuidance(" Available generators : PYTHIA, particleGun");
  fGenCmd->SetParameterName("generator",true);
  fGenCmd->SetDefaultValue("PYTHIA");
  fGenCmd->SetCandidates("PYTHIA particleGun");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01PrimaryGeneratorMessenger::~RE01PrimaryGeneratorMessenger()
{
  delete fGenCmd;
  delete fMydetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01PrimaryGeneratorMessenger
::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fGenCmd )
  { fMyAction->SetHEPEvtGenerator(newValue=="PYTHIA"); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
G4String RE01PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==fGenCmd )
  {
    if(fMyAction->GetHEPEvtGenerator())
    { cv = "PYTHIA"; }
    else
    { cv = "particleGun"; }
  }
  
  return cv;
}

