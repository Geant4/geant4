// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2PhysicsListMessenger.cc,v 1.1 2000-04-17 11:26:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2PhysicsListMessenger.hh"

#include "Em2PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2PhysicsListMessenger::Em2PhysicsListMessenger(Em2PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  gammaCutCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setGCut",this);  
  gammaCutCmd->SetGuidance("Set gamma cut.");
  gammaCutCmd->SetParameterName("Gcut",false);
  gammaCutCmd->SetUnitCategory("Length");
  gammaCutCmd->SetRange("Gcut>0.0");
  gammaCutCmd->AvailableForStates(PreInit,Idle);

  electCutCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setECut",this);  
  electCutCmd->SetGuidance("Set electron cut.");
  electCutCmd->SetParameterName("Ecut",false);
  electCutCmd->SetUnitCategory("Length");
  electCutCmd->SetRange("Ecut>0.0");
  electCutCmd->AvailableForStates(PreInit,Idle);
  
  protoCutCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setPCut",this);  
  protoCutCmd->SetGuidance("Set proton cut.");
  protoCutCmd->SetParameterName("Pcut",false);
  protoCutCmd->SetUnitCategory("Length");
  protoCutCmd->SetRange("Pcut>0.0");
  protoCutCmd->AvailableForStates(PreInit,Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2PhysicsListMessenger::~Em2PhysicsListMessenger()
{
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{       
  if( command == gammaCutCmd )
   { pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == electCutCmd )
   { pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == protoCutCmd )
   { pPhysicsList->SetCutForProton(protoCutCmd->GetNewDoubleValue(newValue));}   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
