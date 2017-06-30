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
// $Id: PrimaryGeneratorMessenger.cc 103795 2017-04-27 13:38:36Z gcosmo $
//
/// \file medical/GammaTherapy/src/PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
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
#include "G4RunManager.hh"
#include "Run.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                           PrimaryGeneratorAction* gen):
  fGen(gen)
{
  fVerbose = gen->GetVerbose();
  if(fVerbose) G4cout << "PrimaryGeneratorMessenger: Construct " << G4endl;

  fBeamXCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamX",this);
  fBeamXCmd->SetGuidance("Set X position of the center of the beam.");
  fBeamXCmd->SetParameterName("beamX",true);
  fBeamXCmd->SetUnitCategory("Length");
  fBeamXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fBeamYCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamY",this);
  fBeamYCmd->SetGuidance("Set Y position of the center of the beam.");
  fBeamYCmd->SetParameterName("beamY",true);
  fBeamYCmd->SetUnitCategory("Length");
  fBeamYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fBeamZCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamZ",this);
  fBeamZCmd->SetGuidance("Set Z of the entry point of the beam.");
  fBeamZCmd->SetParameterName("beamZ",true);
  fBeamZCmd->SetUnitCategory("Length");
  fBeamZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fBeamECmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/beamE",this);
  fBeamECmd->SetGuidance("Set the beam kinetic energy");
  fBeamECmd->SetParameterName("beamE",false);
  fBeamECmd->SetUnitCategory("Energy");
  fBeamECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSigmaXCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaX",this);
  fSigmaXCmd->SetGuidance("Set the beam Gussian width for X");
  fSigmaXCmd->SetParameterName("sigmaX",false);
  fSigmaXCmd->SetUnitCategory("Length");
  fSigmaXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSigmaYCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaY",this);
  fSigmaYCmd->SetGuidance("Set the beam Gussian width for Y");
  fSigmaYCmd->SetParameterName("sigmaY",false);
  fSigmaYCmd->SetUnitCategory("Length");
  fSigmaYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSigmaZCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaZ",this);
  fSigmaZCmd->SetGuidance("Set the beam Gussian width for Y");
  fSigmaZCmd->SetParameterName("sigmaZ",false);
  fSigmaZCmd->SetUnitCategory("Length");
  fSigmaZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSigmaECmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaE",this);
  fSigmaECmd->SetGuidance("Set the beam Gussian width for energy");
  fSigmaECmd->SetParameterName("sigmaE",false);
  fSigmaECmd->SetUnitCategory("Energy");
  fSigmaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fRandCmd = new G4UIcmdWithAString("/testem/gun",this);
  fRandCmd->SetGuidance("Set the name of the random distribution (gauss,flat)");
  fRandCmd->SetParameterName("rand",false);
  fRandCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMaxThetaCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/maxTheta",this);
  fMaxThetaCmd->SetGuidance("Set the beam maxTheta in degrees.");
  fMaxThetaCmd->SetParameterName("maxTheta",false);
  fMaxThetaCmd->SetUnitCategory("Angle");
  fMaxThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fThetaCmd = new G4UIcmdWithADoubleAndUnit("/testem/gun/sigmaTheta",this);
  fThetaCmd->SetGuidance("Set the beam sigmaTheta in degrees.");
  fThetaCmd->SetParameterName("sigmaTheta",false);
  fThetaCmd->SetUnitCategory("Angle");
  fThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fBeamXCmd;
  delete fBeamYCmd;
  delete fBeamZCmd;

  delete fSigmaXCmd;
  delete fSigmaYCmd;
  delete fSigmaZCmd;
  delete fSigmaECmd;

  delete fBeamECmd;
  delete fRandCmd;
  delete fMaxThetaCmd;
  delete fThetaCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if(fVerbose)  
    G4cout << "PrimaryGeneratorMessenger: Next command value = "
           << newValue << G4endl;

  if(command == fBeamXCmd)
    {fGen->SetBeamX(fBeamXCmd->GetNewDoubleValue(newValue));}
  if(command == fBeamYCmd)
    {fGen->SetBeamY(fBeamYCmd->GetNewDoubleValue(newValue));}
  if(command == fBeamZCmd)
    {fGen->SetBeamZ(fBeamZCmd->GetNewDoubleValue(newValue));}
  if(command == fSigmaXCmd)
    {fGen->SetBeamSigmaX(fSigmaXCmd->GetNewDoubleValue(newValue));}
  if(command == fSigmaYCmd)
    {fGen->SetBeamSigmaY(fSigmaYCmd->GetNewDoubleValue(newValue));}
  if(command == fSigmaZCmd)
    {fGen->SetBeamSigmaZ(fSigmaZCmd->GetNewDoubleValue(newValue));}
  if(command == fSigmaECmd)
    {fGen->SetBeamSigmaE(fSigmaECmd->GetNewDoubleValue(newValue));}
  if(command == fBeamECmd) {
    G4double e = fBeamECmd->GetNewDoubleValue(newValue);
    fGen->SetBeamEnergy(e);
  }
  if(command == fMaxThetaCmd)
    {fGen->SetBeamMinCosTheta(
                     std::cos(fMaxThetaCmd->GetNewDoubleValue(newValue)));}
  if(command == fThetaCmd)
    {fGen->SetSigmaTheta(fThetaCmd->GetNewDoubleValue(newValue));}
  if(command == fRandCmd)
    {fGen->SetRandom(newValue);}


  if(fVerbose) G4cout << "PrimaryGeneratorMessenger: O'K " << G4endl;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

