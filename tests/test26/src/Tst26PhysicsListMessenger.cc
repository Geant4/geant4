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
// $Id: Tst26PhysicsListMessenger.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
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
   { pPhysicsList->SetCutForVertex(eCutCmd->GetNewDoubleValue(newValue));}

  if( command == mCutCmd )
   { pPhysicsList->SetCutForMuon(mCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
