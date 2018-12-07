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
/// \file RE05/src/RE05StackingActionMessenger.cc
/// \brief Implementation of the RE05StackingActionMessenger class
//

#include "RE05StackingActionMessenger.hh"
#include "RE05StackingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05StackingActionMessenger::RE05StackingActionMessenger(RE05StackingAction * msa)
: G4UImessenger(),
  fMyAction(msa),
  fMuonCmd(0), fIsoMuonCmd(0), fIsoCmd(0), fRoiCmd(0)
{
  fMuonCmd = new G4UIcmdWithAnInteger("/mydet/reqmuon",this);
  fMuonCmd->SetGuidance("Number of muon for the trigger.");
  fMuonCmd->SetParameterName("N",true);
  fMuonCmd->SetDefaultValue(2);
  fMuonCmd->SetRange("N>=0");

  fIsoMuonCmd = new G4UIcmdWithAnInteger("/mydet/isomuon",this);
  fIsoMuonCmd->SetGuidance("Number of isolated muon for the trigger.");
  fIsoMuonCmd->SetParameterName("N",true);
  fIsoMuonCmd->SetDefaultValue(2);
  fIsoMuonCmd->SetRange("N>=0");

  fIsoCmd = new G4UIcmdWithAnInteger("/mydet/isolation",this);
  fIsoCmd->SetGuidance("Maximum allowed number of hits in tracker");
  fIsoCmd->SetGuidance(" for an isolated muon track (includes hits by muon)");
  fIsoCmd->SetParameterName("N",true);
  fIsoCmd->SetDefaultValue(10);
  fIsoCmd->SetRange("N>=0");

  fRoiCmd = new G4UIcmdWithADoubleAndUnit("/mydet/RoIangle",this);
  fRoiCmd->SetGuidance("Define RoI angle");
  fRoiCmd->SetParameterName("theta",true,true);
  fRoiCmd->SetDefaultUnit("deg");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05StackingActionMessenger::~RE05StackingActionMessenger()
{
  delete fMuonCmd;
  delete fIsoMuonCmd;
  delete fIsoCmd;
  delete fRoiCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05StackingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fMuonCmd )
  { fMyAction->SetNRequestMuon(fMuonCmd->GetNewIntValue(newValue)); }
  else if( command==fIsoMuonCmd )
  { fMyAction->SetNRequestIsoMuon(fIsoMuonCmd->GetNewIntValue(newValue)); }
  else if( command==fIsoCmd )
  { fMyAction->SetNIsolation(fIsoCmd->GetNewIntValue(newValue)); }
  else if( command==fRoiCmd )
  { fMyAction->SetRoIAngle(fRoiCmd->GetNewDoubleValue(newValue)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String RE05StackingActionMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==fMuonCmd )
  { cv = fMuonCmd->ConvertToString(fMyAction->GetNRequestMuon()); }
  else if( command==fIsoMuonCmd )
  { cv = fIsoMuonCmd->ConvertToString(fMyAction->GetNRequestIsoMuon()); }
  else if( command==fIsoCmd )
  { cv = fIsoCmd->ConvertToString(fMyAction->GetNIsolation()); }
  else if( command==fRoiCmd )
  { cv = fRoiCmd->ConvertToString(fMyAction->GetRoIAngle(),"deg"); }
  
  return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
