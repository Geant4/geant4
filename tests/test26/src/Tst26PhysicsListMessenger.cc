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
// $Id: Tst26PhysicsListMessenger.cc,v 1.5 2006-06-29 21:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
// 19-02-03 Rename G4ProductionCuts (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26PhysicsListMessenger.hh"

#include "Tst26PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26PhysicsListMessenger::Tst26PhysicsListMessenger(Tst26PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  wCutCmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/WorldCuts",this);  
  wCutCmd->SetGuidance("Set cuts for the World");
  wCutCmd->SetParameterName("Gcut",false);
  wCutCmd->SetUnitCategory("Length");
  wCutCmd->SetRange("Gcut>0.0");
  wCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  eCutCmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/VertexCuts",this);  
  eCutCmd->SetGuidance("Set cuts for the Vertex Detector");
  eCutCmd->SetParameterName("Ecut",false);
  eCutCmd->SetUnitCategory("Length");
  eCutCmd->SetRange("Ecut>0.0");
  eCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mCutCmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/MuonCuts",this);  
  mCutCmd->SetGuidance("Set cuts for the Muon Detector");
  mCutCmd->SetParameterName("Ecut",false);
  mCutCmd->SetUnitCategory("Length");
  mCutCmd->SetRange("Ecut>0.0");
  mCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  pListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
  pListCmd->SetGuidance("Add modula physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26PhysicsListMessenger::~Tst26PhysicsListMessenger()
{
  delete wCutCmd;
  delete eCutCmd;
  delete mCutCmd;
  delete pListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == wCutCmd )
   { pPhysicsList->SetCutForWorld(wCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == eCutCmd )
   { pPhysicsList->SetCutForVertexDetector(eCutCmd->GetNewDoubleValue(newValue));}

  if( command == mCutCmd )
   { pPhysicsList->SetCutForMuonDetector(mCutCmd->GetNewDoubleValue(newValue));}

  if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
