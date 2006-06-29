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

//---------------------------------------------------------------------------
//
// ClassName:   PrimaryGeneratorMessenger
//
// Description: Definition of physics list parameters
//
// Author:      V.Ivanchenko 26/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "Histo.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                PrimaryGeneratorAction* gen):
  theGen(gen)
{
  G4cout << "PrimaryGeneratorMessenger: Construct " << G4endl;

  beamXCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamX",this);
  beamXCmd->SetGuidance("Set X position of the center of the beam.");
  beamXCmd->SetParameterName("beamX",true);
  beamXCmd->SetUnitCategory("Length");
  beamXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamYCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamY",this);
  beamYCmd->SetGuidance("Set Y position of the center of the beam.");
  beamYCmd->SetParameterName("beamY",true);
  beamYCmd->SetUnitCategory("Length");
  beamYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamZCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamZ",this);
  beamZCmd->SetGuidance("Set Z of the entry point of the beam.");
  beamZCmd->SetParameterName("beamZ",true);
  beamZCmd->SetUnitCategory("Length");
  beamZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sigmaXCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaX",this);
  sigmaXCmd->SetGuidance("Set the beam Gussian width for X");
  sigmaXCmd->SetParameterName("sigmaX",false);
  sigmaXCmd->SetUnitCategory("Length");
  sigmaXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sigmaYCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaY",this);
  sigmaYCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaYCmd->SetParameterName("sigmaY",false);
  sigmaYCmd->SetUnitCategory("Length");
  sigmaYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sigmaZCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaZ",this);
  sigmaZCmd->SetGuidance("Set the beam Gussian width for Y");
  sigmaZCmd->SetParameterName("sigmaZ",false);
  sigmaZCmd->SetUnitCategory("Length");
  sigmaZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sigmaECmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaE",this);
  sigmaECmd->SetGuidance("Set the beam Gussian width for energy");
  sigmaECmd->SetParameterName("sigmaE",false);
  sigmaECmd->SetUnitCategory("Energy");
  sigmaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamECmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamE",this);
  beamECmd->SetGuidance("Set the beam kinetic energy");
  beamECmd->SetParameterName("beamE",false);
  beamECmd->SetUnitCategory("Energy");
  beamECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  randCmd = new G4UIcmdWithAString("/testem/gun",this);
  randCmd->SetGuidance("Set the name of the random distribution (gauss,flat)");
  randCmd->SetParameterName("rand",false);
  randCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  maxThetaCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/maxTheta",this);
  maxThetaCmd->SetGuidance("Set the beam maxTheta in degrees.");
  maxThetaCmd->SetParameterName("maxTheta",false);
  maxThetaCmd->SetUnitCategory("Angle");
  maxThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  sThetaCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaTheta",this);
  sThetaCmd->SetGuidance("Set the beam sigmaTheta in degrees.");
  sThetaCmd->SetParameterName("sigmaTheta",false);
  sThetaCmd->SetUnitCategory("Angle");
  sThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
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
  delete sThetaCmd;
  delete randCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                                 G4String newValue)
{

  if(1 < (Histo::GetPointer())->GetVerbose()) {
    G4cout << "PrimaryGeneratorMessenger: Next command value = "
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
    Histo* theHisto = Histo::GetPointer();
    if(theHisto->GetMaxEnergy() == 0.0) theHisto->SetMaxEnergy(e);
  }
  if(command == maxThetaCmd)
    {theGen->SetBeamMinCosTheta(std::cos(maxThetaCmd->GetNewDoubleValue(newValue)));}
  if(command == sThetaCmd)
    {theGen->SetSigmaTheta(sThetaCmd->GetNewDoubleValue(newValue));}
  if(command == randCmd)
    {theGen->SetRandom(newValue);}


  if(1 < (Histo::GetPointer())->GetVerbose())
    {G4cout << "PrimaryGeneratorMessenger: O'K " << G4endl;}
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

