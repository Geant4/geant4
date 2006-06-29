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
// $Id: A01MagneticFieldMessenger.cc,v 1.4 2006-06-29 16:32:59 gunter Exp $
// --------------------------------------------------------------
//
#include "A01MagneticFieldMessenger.hh"
#include "A01MagneticField.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

A01MagneticFieldMessenger::A01MagneticFieldMessenger(A01MagneticField * mpga)
:target(mpga)
{
  fieldCmd = new G4UIcmdWithADoubleAndUnit("/mydet/fieldValue",this);
  fieldCmd->SetGuidance("Field strength");
  fieldCmd->SetParameterName("field",true);
  fieldCmd->SetDefaultValue(1.);
  fieldCmd->SetDefaultUnit("tesla");
}

A01MagneticFieldMessenger::~A01MagneticFieldMessenger()
{
  delete fieldCmd;
}

void A01MagneticFieldMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fieldCmd )
  { target->SetField(fieldCmd->GetNewDoubleValue(newValue)); }
}

G4String A01MagneticFieldMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==fieldCmd )
  { cv = fieldCmd->ConvertToString(target->GetField(),"tesla"); }

  return cv;
}

