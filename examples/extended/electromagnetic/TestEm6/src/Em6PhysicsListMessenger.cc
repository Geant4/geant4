// Em6PhysicsListMessenger.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6PhysicsListMessenger.hh"

#include "Em6PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6PhysicsListMessenger::Em6PhysicsListMessenger(Em6PhysicsList* pPhys)
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

  GammaToMuPairFac=new G4UIcmdWithADouble("/run/process/setGammaToMuPairFac",this);
  GammaToMuPairFac->SetGuidance(
         "Set factor to artificially increase the GammaToMuPair cross section");
  GammaToMuPairFac->SetParameterName("GammaToMuPairFac",false);
  GammaToMuPairFac->SetRange("GammaToMuPairFac>0.0");
  GammaToMuPairFac->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6PhysicsListMessenger::~Em6PhysicsListMessenger()
{
  delete gammaCutCmd;
  delete electCutCmd;
  delete protoCutCmd;
  delete GammaToMuPairFac;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == gammaCutCmd )
   { pPhysicsList->SetCutForGamma(gammaCutCmd->GetNewDoubleValue(newValue));}

  if( command == electCutCmd )
   { pPhysicsList->SetCutForElectron(electCutCmd->GetNewDoubleValue(newValue));}

  if( command == protoCutCmd )
   { pPhysicsList->SetCutForProton(protoCutCmd->GetNewDoubleValue(newValue));}

  if( command == GammaToMuPairFac )
   { pPhysicsList->setGammaToMuPairFac(
                              GammaToMuPairFac->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
