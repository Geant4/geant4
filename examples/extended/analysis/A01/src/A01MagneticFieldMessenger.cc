// $Id: A01MagneticFieldMessenger.cc,v 1.1 2002-11-13 07:23:50 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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

