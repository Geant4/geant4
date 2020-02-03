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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorMessenger  ------
//           by G.Santin, F.Longo & R.Giannitrapani (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelPrimaryGeneratorMessenger.hh"

#include "GammaRayTelPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorMessenger::GammaRayTelPrimaryGeneratorMessenger
(GammaRayTelPrimaryGeneratorAction* GammaRayTelGun)
  :GammaRayTelAction(GammaRayTelGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SourceTypeCmd = new G4UIcmdWithAnInteger("/gun/sourceType",this);
  SourceTypeCmd->SetGuidance("Select the type of incident flux.");
  SourceTypeCmd->SetGuidance("  Choice : 0(default), 1(isotropic), 2(wide parallel beam)");
  SourceTypeCmd->SetParameterName("choice",true);
  SourceTypeCmd->SetDefaultValue((G4int)0);
  SourceTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  VertexRadiusCmd = new G4UIcmdWithADoubleAndUnit("/gun/vertexRadius",this);
  VertexRadiusCmd->SetGuidance("Radius (and unit) of sphere for vertices of incident flux.");
  VertexRadiusCmd->SetParameterName("choice",true);
  VertexRadiusCmd->SetDefaultValue((G4double)1.*cm);
  VertexRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SpectrumTypeCmd = new G4UIcmdWithAnInteger("/gun/spectrumType",this);
  SpectrumTypeCmd->SetGuidance("Select the type of incident spectrum.");
  SpectrumTypeCmd->SetGuidance("  Choice : 0(default), 1(), 2(E^{-gamma}), 3()");
  SpectrumTypeCmd->SetParameterName("choice",true);
  SpectrumTypeCmd->SetDefaultValue((G4int)0);
  SpectrumTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  SourceGenCmd = new G4UIcmdWithABool("/gun/sourceGen",this);
  SourceGenCmd->SetGuidance("Select the native Generation");
  SourceGenCmd->SetGuidance("  Choice : true(native), false(GPS)");
  SourceGenCmd->SetParameterName("choice",true);
  SourceGenCmd->SetDefaultValue((G4bool)true);
  SourceGenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorMessenger::~GammaRayTelPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete SourceTypeCmd;
  delete VertexRadiusCmd;
  delete SpectrumTypeCmd;
  delete SourceGenCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd )
    { GammaRayTelAction->SetRndmFlag(newValue);}
  
  if( command == SourceTypeCmd )
    { GammaRayTelAction->SetSourceType(SourceTypeCmd->GetNewIntValue(newValue));}
  
  if( command == VertexRadiusCmd )
    { GammaRayTelAction->SetVertexRadius(VertexRadiusCmd->GetNewDoubleValue(newValue));}
  
  if( command == SpectrumTypeCmd )
    { GammaRayTelAction->SetSpectrumType(SpectrumTypeCmd->GetNewIntValue(newValue));}

  if( command == SourceGenCmd )
    { GammaRayTelAction->SetSourceGen(SourceGenCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












