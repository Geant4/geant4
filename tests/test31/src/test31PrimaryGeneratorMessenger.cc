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
#define test31PrimaryGeneratorMessenger_CPP 

//---------------------------------------------------------------------------
//
// ClassName:   test31PrimaryGeneratorMessenger
//  
// Description: Definition of physics list parameters
//
// Author:      V.Ivanchenko 26/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31PrimaryGeneratorMessenger.hh"
#include "test31PrimaryGeneratorAction.hh"
#include "test31Histo.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PrimaryGeneratorMessenger::test31PrimaryGeneratorMessenger(
                                test31PrimaryGeneratorAction* gen):
  theGen(gen)
{
  G4cout << "test31PrimaryGeneratorMessenger: Construct " << G4endl;

  beamXCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/beamX",this);
  beamXCmd->SetGuidance("Set X position of the center of the beam.");
  beamXCmd->SetParameterName("beamX",true);
  beamXCmd->SetUnitCategory("Length");
  beamXCmd->AvailableForStates(PreInit,Idle);

  beamYCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/beamY",this);
  beamYCmd->SetGuidance("Set Y position of the center of the beam.");
  beamYCmd->SetParameterName("beamY",true);
  beamYCmd->SetUnitCategory("Length");
  beamYCmd->AvailableForStates(PreInit,Idle);

  beamZCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/beamZ",this);
  beamZCmd->SetGuidance("Set Z of the entry point of the beam.");
  beamZCmd->SetParameterName("beamZ",true);
  beamZCmd->SetUnitCategory("Length");
  beamZCmd->AvailableForStates(PreInit,Idle);

  sigmaXCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/sigmaX",this);
  sigmaXCmd->SetGuidance("Set the beam Gussian width for X");
  sigmaXCmd->SetParameterName("sigmaX",false);
  sigmaXCmd->SetUnitCategory("Length");
  sigmaXCmd->AvailableForStates(PreInit,Idle);

  sigmaYCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/sigmaY",this);
  sigmaYCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaYCmd->SetParameterName("sigmaY",false);
  sigmaYCmd->SetUnitCategory("Length");
  sigmaYCmd->AvailableForStates(PreInit,Idle);

  sigmaZCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/sigmaZ",this);
  sigmaZCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaZCmd->SetParameterName("sigmaZ",false);
  sigmaZCmd->SetUnitCategory("Length");
  sigmaZCmd->AvailableForStates(PreInit,Idle);

  sigmaECmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/sigmaE",this);
  sigmaECmd->SetGuidance("Set the beam Gussian width for energy");
  sigmaECmd->SetParameterName("sigmaE",false);
  sigmaECmd->SetUnitCategory("Energy");
  sigmaECmd->AvailableForStates(PreInit,Idle);

  beamECmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/beamE",this);
  beamECmd->SetGuidance("Set the beam kinetic energy");
  beamECmd->SetParameterName("beamE",false);
  beamECmd->SetUnitCategory("Energy");
  beamECmd->AvailableForStates(PreInit,Idle);

  beamBetaCmd = new G4UIcmdWithADouble("/test31/gun/beamBeta",this);
  beamBetaCmd->SetGuidance("Set the beam velocity");
  beamBetaCmd->SetParameterName("beamBeta",false);
  beamBetaCmd->AvailableForStates(PreInit,Idle);

  sigmaBetaCmd = new G4UIcmdWithADouble("/test31/gun/sigmaBeta",this);
  sigmaBetaCmd->SetGuidance("Set the sigma velocity");
  sigmaBetaCmd->SetParameterName("sigmaBeta",false);
  sigmaBetaCmd->AvailableForStates(PreInit,Idle);

  randCmd = new G4UIcmdWithAString("/test31/gun/random",this);
  randCmd->SetGuidance("Set the name of the random distribution (gauss,flatE,flatBeta)");
  randCmd->SetParameterName("rand",false);
  randCmd->AvailableForStates(PreInit,Idle);

  partCmd = new G4UIcmdWithAString("/test31/gun/particle",this);
  partCmd->SetGuidance("Set the name of the particle");
  partCmd->SetParameterName("part",false);
  partCmd->AvailableForStates(PreInit,Idle);

  maxThetaCmd = new G4UIcmdWithADoubleAndUnit("/test31/gun/maxTheta",this);
  maxThetaCmd->SetGuidance("Set the beam maxTheta in degrees.");
  maxThetaCmd->SetParameterName("maxTheta",false);
  maxThetaCmd->SetUnitCategory("Angle");
  maxThetaCmd->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PrimaryGeneratorMessenger::~test31PrimaryGeneratorMessenger()
{
  delete beamXCmd;
  delete beamYCmd;
  delete beamZCmd;
  delete sigmaXCmd;
  delete sigmaYCmd;
  delete sigmaZCmd;
  delete sigmaECmd;
  delete beamECmd;
  delete maxThetaCmd;
  delete partCmd;
  delete beamBetaCmd;
  delete sigmaBetaCmd;
  delete randCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void test31PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                                 G4String newValue)
{
 
  if(1 < theGen->GetVerbose()) {
    G4cout << "test31PrimaryGeneratorMessenger: Next command value = " 
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
  if(command == beamECmd) {
    G4double e = beamECmd->GetNewDoubleValue(newValue);
    theGen->SetBeamEnergy(e);
    test31Histo* theHisto = test31Histo::GetPointer();
    if(theHisto->GetMaxEnergy() == 0.0) theHisto->SetMaxEnergy(e);
  }
  if(command == maxThetaCmd)
    {theGen->SetBeamMinCosTheta(cos(maxThetaCmd->GetNewDoubleValue(newValue)));}
  if(command == partCmd)
    {(G4UImanager::GetUIpointer())->ApplyCommand("/gun/particle "+newValue);}
  if(command == beamBetaCmd)
    {theGen->SetBeamBeta(beamBetaCmd->GetNewDoubleValue(newValue));}
  if(command == sigmaBetaCmd)
    {theGen->SetSigmaBeta(sigmaBetaCmd->GetNewDoubleValue(newValue));}
  if(command == randCmd)
    {theGen->SetRandom(newValue);}


  if(0 < theGen->GetVerbose())
    {G4cout << "test31PrimaryGeneratorMessenger: O'K " << G4endl;}
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

