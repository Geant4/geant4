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
// $Id: NTSTPhysicsListMessenger.cc,v 1.2 2003-12-09 15:35:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NTSTPhysicsListMessenger.hh"

#include "NTSTPhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTPhysicsListMessenger::NTSTPhysicsListMessenger(NTSTPhysicsList * List)
:NTSTList(List)
{
  ProcessCmd = new G4UIcmdWithoutParameter("/run/setEmProcess",this);
  ProcessCmd->SetGuidance("select electromagnetic processes");
  ProcessCmd->AvailableForStates(G4State_PreInit);  // New-state
  // ProcessCmd->AvailableForStates(PreInit);  // Old-state

  useBgsTranCmd = new G4UIcmdWithABool("/run/useBgsTran",this);
  useBgsTranCmd->SetGuidance( "True if BgsTransportation to be used" );
  useBgsTranCmd->AvailableForStates(G4State_PreInit);  // New-state
  // useBgsTranCmd->AvailableForStates(PreInit);  // Old-state

  MinimumEnergyCutCmd = new G4UIcmdWithADoubleAndUnit("/run/minEcut", this);
  MinimumEnergyCutCmd->SetGuidance("ParticleWithCuts minimum energy cut (with unit)" );
  MinimumEnergyCutCmd->AvailableForStates(G4State_PreInit);  // New-state
  // MinimumEnergyCutCmd->AvailableForStates(PreInit);  // Old-state

  MaximumEnergyCutCmd = new G4UIcmdWithADoubleAndUnit("/run/maxEcut", this);
  MaximumEnergyCutCmd->SetGuidance("ParticleWithCuts maximum energy cut (with unit)" );
  MaximumEnergyCutCmd->AvailableForStates(G4State_PreInit);  // New-state
  // MaximumEnergyCutCmd->AvailableForStates(PreInit);  // Old-state

  CutCmd = new G4UIcmdWithADoubleAndUnit("/run/cut", this);
  CutCmd->SetGuidance("new cut value (with unit)" );
  CutCmd->AvailableForStates(G4State_PreInit);  // New-state
  // CutCmd->AvailableForStates(PreInit);  // Old-state

  LooperCutCmd  = new G4UIcmdWithADoubleAndUnit("/run/looperCut", this);
  LooperCutCmd->SetGuidance("Kill loopers below this cut (with unit)" );
  LooperCutCmd->AvailableForStates(G4State_PreInit);  // New-state
  // LooperCutCmd->AvailableForStates(PreInit);  // Old-state
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTPhysicsListMessenger::~NTSTPhysicsListMessenger()
{
  delete ProcessCmd;
  delete useBgsTranCmd;
  delete MinimumEnergyCutCmd;
  delete MaximumEnergyCutCmd;
  delete CutCmd;
  delete LooperCutCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void NTSTPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String s)
{
  if(command == ProcessCmd) {
    NTSTList->SetStatusEmProcess();
  } else if (command == useBgsTranCmd) {
    NTSTList->SetBgsTran( useBgsTranCmd->GetNewBoolValue(s) );
  } else if (command == MinimumEnergyCutCmd) {
    NTSTList->SetMinimumEnergyCut( MinimumEnergyCutCmd->GetNewDoubleValue(s));
  } else if (command == MaximumEnergyCutCmd) {
    NTSTList->SetMaximumEnergyCut( MaximumEnergyCutCmd->GetNewDoubleValue(s));
  } else if (command == CutCmd) {
    NTSTList->SetLengthCut( CutCmd->GetNewDoubleValue(s) );
  } else if (command == LooperCutCmd) {
    NTSTList->SetLooperCut( LooperCutCmd->GetNewDoubleValue(s) );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







