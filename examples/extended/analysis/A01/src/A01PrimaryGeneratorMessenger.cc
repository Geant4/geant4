// $Id: A01PrimaryGeneratorMessenger.cc,v 1.1 2002-11-13 07:24:05 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "A01PrimaryGeneratorMessenger.hh"
#include "A01PrimaryGeneratorAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"

A01PrimaryGeneratorMessenger::A01PrimaryGeneratorMessenger(A01PrimaryGeneratorAction * mpga)
:target(mpga)
{
  momentumCmd = new G4UIcmdWithADoubleAndUnit("/mydet/momentum",this);
  momentumCmd->SetGuidance("Mean momentum of primaries");
  momentumCmd->SetParameterName("p",true);
  momentumCmd->SetRange("p>=0.");
  momentumCmd->SetDefaultValue(1.);
  momentumCmd->SetDefaultUnit("GeV");

  sigmaMomCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sigmaMomentum",this);
  sigmaMomCmd->SetGuidance("Sigma momentum of primaries");
  sigmaMomCmd->SetParameterName("p",true);
  sigmaMomCmd->SetRange("p>=0.");
  sigmaMomCmd->SetDefaultValue(1.);
  sigmaMomCmd->SetDefaultUnit("GeV");

  sigmaAngCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sigmaAngle",this);
  sigmaAngCmd->SetGuidance("sigma angle divergence of primaries");
  sigmaAngCmd->SetParameterName("t",true);
  sigmaAngCmd->SetRange("t>=0.");
  sigmaAngCmd->SetDefaultValue(1.);
  sigmaAngCmd->SetDefaultUnit("deg");

  randomCmd = new G4UIcmdWithABool("/mydet/randomizePrimary",this);
  randomCmd->SetGuidance("Boolean flag for randomizing primary particle types.");
  randomCmd->SetGuidance("In case this flag is false, you can select the primary particle");
  randomCmd->SetGuidance("  with /gun/particle command.");
  randomCmd->SetParameterName("flg",true);
  randomCmd->SetDefaultValue(true);
}

A01PrimaryGeneratorMessenger::~A01PrimaryGeneratorMessenger()
{
  delete momentumCmd;
  delete sigmaMomCmd;
  delete sigmaAngCmd;
  delete randomCmd;
}

void A01PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==momentumCmd )
  { target->SetMomentum(momentumCmd->GetNewDoubleValue(newValue)); }
  if( command==sigmaMomCmd )
  { target->SetSigmaMomentum(sigmaMomCmd->GetNewDoubleValue(newValue)); }
  if( command==sigmaAngCmd )
  { target->SetSigmaAngle(sigmaAngCmd->GetNewDoubleValue(newValue)); }
  if( command==randomCmd )
  { target->SetRandomize(randomCmd->GetNewBoolValue(newValue)); }
}

G4String A01PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==momentumCmd )
  { cv = momentumCmd->ConvertToString(target->GetMomentum(),"GeV"); }
  if( command==sigmaMomCmd )
  { cv = sigmaMomCmd->ConvertToString(target->GetSigmaMomentum(),"GeV"); }
  if( command==sigmaAngCmd )
  { cv = sigmaAngCmd->ConvertToString(target->GetSigmaAngle(),"deg"); }
  if( command==randomCmd )
  { cv = randomCmd->ConvertToString(target->GetRandomize()); }

  return cv;
}

