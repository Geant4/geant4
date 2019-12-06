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
/// \file ExN04PrimaryGeneratorMessenger.cc
/// \brief Implementation of the ExN04PrimaryGeneratorMessenger class
//

#include "ExN04PrimaryGeneratorMessenger.hh"
#include "ExN04PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

ExN04PrimaryGeneratorMessenger::ExN04PrimaryGeneratorMessenger(ExN04PrimaryGeneratorAction * mpga)
:myAction(mpga)
{
  mydetDirectory = new G4UIdirectory("/mydet/");
  mydetDirectory->SetGuidance("ExN04 detector control commands.");

  genCmd = new G4UIcmdWithAString("/mydet/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  genCmd->SetGuidance(" Available generators : PYTHIA, particleGun");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("PYTHIA");
  genCmd->SetCandidates("PYTHIA particleGun");
}

ExN04PrimaryGeneratorMessenger::~ExN04PrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
}

void ExN04PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  { myAction->SetHEPEvtGenerator(newValue=="PYTHIA"); }
}

G4String ExN04PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->GetHEPEvtGenerator())
    { cv = "PYTHIA"; }
    else
    { cv = "particleGun"; }
  }
  
  return cv;
}

