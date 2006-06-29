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
// $Id: RE01PrimaryGeneratorMessenger.cc,v 1.2 2006-06-29 17:44:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#include "RE01PrimaryGeneratorMessenger.hh"
#include "RE01PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

RE01PrimaryGeneratorMessenger::RE01PrimaryGeneratorMessenger(RE01PrimaryGeneratorAction * mpga)
:myAction(mpga)
{
  mydetDirectory = new G4UIdirectory("/mydet/");
  mydetDirectory->SetGuidance("RE01 detector control commands.");

  genCmd = new G4UIcmdWithAString("/mydet/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  genCmd->SetGuidance(" Available generators : PYTHIA, particleGun");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("PYTHIA");
  genCmd->SetCandidates("PYTHIA particleGun");
}

RE01PrimaryGeneratorMessenger::~RE01PrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
}

void RE01PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  { myAction->SetHEPEvtGenerator(newValue=="PYTHIA"); }
}

G4String RE01PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
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

