#define hTestPrimaryGeneratorMessenger_CPP 

//---------------------------------------------------------------------------
//
// ClassName:   hTestPrimaryGeneratorMessenger
//  
// Description: Definition of physics list parameters
//
// Author:      V.Ivanchenko 26/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestPrimaryGeneratorMessenger.hh"
#include "hTestPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorMessenger::hTestPrimaryGeneratorMessenger(
                                hTestPrimaryGeneratorAction* gen):
  theGen(gen)
{
  G4cout << "hTestPrimaryGeneratorMessenger: Construct " << G4endl;

  beamXCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/beamX",this);
  beamXCmd->SetGuidance("Set X position of the center of the beam.");
  beamXCmd->SetParameterName("beamX",true);
  beamXCmd->SetUnitCategory("Length");
  beamXCmd->AvailableForStates(PreInit,Idle);

  beamYCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/beamY",this);
  beamYCmd->SetGuidance("Set Y position of the center of the beam.");
  beamYCmd->SetParameterName("beamY",true);
  beamYCmd->SetUnitCategory("Length");
  beamYCmd->AvailableForStates(PreInit,Idle);

  beamZCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/beamZ",this);
  beamZCmd->SetGuidance("Set Z of the entry point of the beam.");
  beamZCmd->SetParameterName("beamZ",true);
  beamZCmd->SetUnitCategory("Length");
  beamZCmd->AvailableForStates(PreInit,Idle);

  sigmaXCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/sigmaX",this);
  sigmaXCmd->SetGuidance("Set the beam Gussian width for X");
  sigmaXCmd->SetParameterName("sigmaX",false);
  sigmaXCmd->SetUnitCategory("Length");
  sigmaXCmd->AvailableForStates(PreInit,Idle);

  sigmaYCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/sigmaY",this);
  sigmaYCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaYCmd->SetParameterName("sigmaY",false);
  sigmaYCmd->SetUnitCategory("Length");
  sigmaYCmd->AvailableForStates(PreInit,Idle);

  sigmaZCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/sigmaZ",this);
  sigmaZCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaZCmd->SetParameterName("sigmaZ",false);
  sigmaZCmd->SetUnitCategory("Length");
  sigmaZCmd->AvailableForStates(PreInit,Idle);

  sigmaECmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/sigmaE",this);
  sigmaECmd->SetGuidance("Set the beam Gussian width for energy");
  sigmaECmd->SetParameterName("sigmaE",false);
  sigmaECmd->SetUnitCategory("Energy");
  sigmaECmd->AvailableForStates(PreInit,Idle);

  maxThetaCmd = new G4UIcmdWithADoubleAndUnit("/hTest/gun/maxTheta",this);
  maxThetaCmd->SetGuidance("Set the beam maxTheta in degrees.");
  maxThetaCmd->SetParameterName("maxTheta",false);
  maxThetaCmd->SetUnitCategory("Angle");
  maxThetaCmd->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorMessenger::~hTestPrimaryGeneratorMessenger()
{
  delete beamXCmd;
  delete beamYCmd;
  delete beamZCmd;
  delete sigmaXCmd;
  delete sigmaYCmd;
  delete sigmaZCmd;
  delete sigmaECmd;
  delete maxThetaCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void hTestPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                                 G4String newValue)
{
 
  if(1 < theGen->GetVerbose()) {
    G4cout << "hTestPrimaryGeneratorMessenger: Next command value = " 
           << newValue << G4endl;
  }

  if(command == beamXCmd)
    {theGen->SetBeamX(beamXCmd->GetNewDoubleValue(newValue));}
  if(command == beamYCmd)
    {theGen->SetBeamY(beamYCmd->GetNewDoubleValue(newValue));}
  if(command == beamZCmd) 
    {theGen->SetBeamZ(beamZCmd->GetNewDoubleValue(newValue));}
  if(command == sigmaXCmd)
    {theGen->SetBeamSigmaX(sigmaXCmd->GetNewDoubleValue(newValue));}
  if(command == sigmaYCmd)
    {theGen->SetBeamSigmaY(sigmaYCmd->GetNewDoubleValue(newValue));}
  if(command == sigmaZCmd)
    {theGen->SetBeamSigmaZ(sigmaZCmd->GetNewDoubleValue(newValue));}
  if(command == sigmaECmd)
    {theGen->SetBeamSigmaE(sigmaECmd->GetNewDoubleValue(newValue));}
  if(command == maxThetaCmd)
    {theGen->SetBeamMinCosTheta(cos(maxThetaCmd->GetNewDoubleValue(newValue)));}
  if(0 < theGen->GetVerbose())
    {G4cout << "hTestPrimaryGeneratorMessenger: O'K " << G4endl;}
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

