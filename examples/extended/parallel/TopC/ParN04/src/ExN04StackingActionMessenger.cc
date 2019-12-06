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
/// \file ExN04StackingActionMessenger.cc
/// \brief Implementation of the ExN04StackingActionMessenger class
//

#include "ExN04StackingActionMessenger.hh"
#include "ExN04StackingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

ExN04StackingActionMessenger::ExN04StackingActionMessenger(ExN04StackingAction * msa)
:myAction(msa)
{
  muonCmd = new G4UIcmdWithAnInteger("/mydet/reqmuon",this);
  muonCmd->SetGuidance("Number of muon for the trigger.");
  muonCmd->SetParameterName("N",true);
  muonCmd->SetDefaultValue(2);
  muonCmd->SetRange("N>=0");

  isomuonCmd = new G4UIcmdWithAnInteger("/mydet/isomuon",this);
  isomuonCmd->SetGuidance("Number of isolated muon for the trigger.");
  isomuonCmd->SetParameterName("N",true);
  isomuonCmd->SetDefaultValue(2);
  isomuonCmd->SetRange("N>=0");

  isoCmd = new G4UIcmdWithAnInteger("/mydet/isolation",this);
  isoCmd->SetGuidance("Maximum allowed number of hits in tracker");
  isoCmd->SetGuidance(" for an isolated muon track (includes hits by muon)");
  isoCmd->SetParameterName("N",true);
  isoCmd->SetDefaultValue(10);
  isoCmd->SetRange("N>=0");

  roiCmd = new G4UIcmdWithADoubleAndUnit("/mydet/RoIangle",this);
  roiCmd->SetGuidance("Define RoI angle");
  roiCmd->SetParameterName("theta",true,true);
  roiCmd->SetDefaultUnit("deg");
}

ExN04StackingActionMessenger::~ExN04StackingActionMessenger()
{
  delete muonCmd;
  delete isomuonCmd;
  delete isoCmd;
  delete roiCmd;
}

void ExN04StackingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==muonCmd )
  { myAction->SetNRequestMuon(muonCmd->GetNewIntValue(newValue)); }
  else if( command==isomuonCmd )
  { myAction->SetNRequestIsoMuon(isomuonCmd->GetNewIntValue(newValue)); }
  else if( command==isoCmd )
  { myAction->SetNIsolation(isoCmd->GetNewIntValue(newValue)); }
  else if( command==roiCmd )
  { myAction->SetRoIAngle(roiCmd->GetNewDoubleValue(newValue)); }
}

G4String ExN04StackingActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==muonCmd )
  { cv = muonCmd->ConvertToString(myAction->GetNRequestMuon()); }
  else if( command==isomuonCmd )
  { cv = isomuonCmd->ConvertToString(myAction->GetNRequestIsoMuon()); }
  else if( command==isoCmd )
  { cv = isoCmd->ConvertToString(myAction->GetNIsolation()); }
  else if( command==roiCmd )
  { cv = roiCmd->ConvertToString(myAction->GetRoIAngle(),"deg"); }
  
  return cv;
}

