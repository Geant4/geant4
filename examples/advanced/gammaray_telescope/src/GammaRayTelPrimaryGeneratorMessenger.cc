//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: GammaRayTelPrimaryGeneratorMessenger.cc,v 1.3 2001-07-11 09:56:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

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
  RndmCmd->AvailableForStates(PreInit,Idle);
  
  SourceTypeCmd = new G4UIcmdWithAnInteger("/gun/sourceType",this);
  SourceTypeCmd->SetGuidance("Select the type of incident flux.");
  SourceTypeCmd->SetGuidance("  Choice : 0(default), 1(isotropic), 2(wide parallel beam)");
  SourceTypeCmd->SetParameterName("choice",true);
  SourceTypeCmd->SetDefaultValue((G4int)0);
  SourceTypeCmd->AvailableForStates(PreInit,Idle);

  VertexRadiusCmd = new G4UIcmdWithADoubleAndUnit("/gun/vertexRadius",this);
  VertexRadiusCmd->SetGuidance("Radius (and unit) of sphere for vertices of incident flux.");
  VertexRadiusCmd->SetParameterName("choice",true);
  VertexRadiusCmd->SetDefaultValue((G4double)1.*cm);
  VertexRadiusCmd->AvailableForStates(PreInit,Idle);
  
  SpectrumTypeCmd = new G4UIcmdWithAnInteger("/gun/spectrumType",this);
  SpectrumTypeCmd->SetGuidance("Select the type of incident spectrum.");
  SpectrumTypeCmd->SetGuidance("  Choice : 0(default), 1(), 2(E^{-gamma}), 3()");
  SpectrumTypeCmd->SetParameterName("choice",true);
  SpectrumTypeCmd->SetDefaultValue((G4int)0);
  SpectrumTypeCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorMessenger::~GammaRayTelPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete SourceTypeCmd;
  delete VertexRadiusCmd;
  delete SpectrumTypeCmd;
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












