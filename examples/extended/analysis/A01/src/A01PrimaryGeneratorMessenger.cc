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
/// \file analysis/A01/src/A01PrimaryGeneratorMessenger.cc
/// \brief Implementation of the A01PrimaryGeneratorMessenger class
//
// $Id$
// --------------------------------------------------------------
//
#include "A01PrimaryGeneratorMessenger.hh"
#include "A01PrimaryGeneratorAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"

A01PrimaryGeneratorMessenger::A01PrimaryGeneratorMessenger(A01PrimaryGeneratorAction * mpga)
:fTarget(mpga)
{
  fMomentumCmd = new G4UIcmdWithADoubleAndUnit("/mydet/momentum",this);
  fMomentumCmd->SetGuidance("Mean momentum of primaries");
  fMomentumCmd->SetParameterName("p",true);
  fMomentumCmd->SetRange("p>=0.");
  fMomentumCmd->SetDefaultValue(1.);
  fMomentumCmd->SetDefaultUnit("GeV");

  fSigmaMomCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sigmaMomentum",this);
  fSigmaMomCmd->SetGuidance("Sigma momentum of primaries");
  fSigmaMomCmd->SetParameterName("p",true);
  fSigmaMomCmd->SetRange("p>=0.");
  fSigmaMomCmd->SetDefaultValue(1.);
  fSigmaMomCmd->SetDefaultUnit("GeV");

  fSigmaAngCmd = new G4UIcmdWithADoubleAndUnit("/mydet/sigmaAngle",this);
  fSigmaAngCmd->SetGuidance("sigma angle divergence of primaries");
  fSigmaAngCmd->SetParameterName("t",true);
  fSigmaAngCmd->SetRange("t>=0.");
  fSigmaAngCmd->SetDefaultValue(1.);
  fSigmaAngCmd->SetDefaultUnit("deg");

  fRandomCmd = new G4UIcmdWithABool("/mydet/randomizePrimary",this);
  fRandomCmd->SetGuidance("Boolean flag for randomizing primary particle types.");
  fRandomCmd->SetGuidance("In case this flag is false, you can select the primary particle");
  fRandomCmd->SetGuidance("  with /gun/particle command.");
  fRandomCmd->SetParameterName("flg",true);
  fRandomCmd->SetDefaultValue(true);
}

A01PrimaryGeneratorMessenger::~A01PrimaryGeneratorMessenger()
{
  delete fMomentumCmd;
  delete fSigmaMomCmd;
  delete fSigmaAngCmd;
  delete fRandomCmd;
}

void A01PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fMomentumCmd )
  { fTarget->SetMomentum(fMomentumCmd->GetNewDoubleValue(newValue)); }
  if( command==fSigmaMomCmd )
  { fTarget->SetSigmaMomentum(fSigmaMomCmd->GetNewDoubleValue(newValue)); }
  if( command==fSigmaAngCmd )
  { fTarget->SetSigmaAngle(fSigmaAngCmd->GetNewDoubleValue(newValue)); }
  if( command==fRandomCmd )
  { fTarget->SetRandomize(fRandomCmd->GetNewBoolValue(newValue)); }
}

G4String A01PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==fMomentumCmd )
  { cv = fMomentumCmd->ConvertToString(fTarget->GetMomentum(),"GeV"); }
  if( command==fSigmaMomCmd )
  { cv = fSigmaMomCmd->ConvertToString(fTarget->GetSigmaMomentum(),"GeV"); }
  if( command==fSigmaAngCmd )
  { cv = fSigmaAngCmd->ConvertToString(fTarget->GetSigmaAngle(),"deg"); }
  if( command==fRandomCmd )
  { cv = fRandomCmd->ConvertToString(fTarget->GetRandomize()); }

  return cv;
}

