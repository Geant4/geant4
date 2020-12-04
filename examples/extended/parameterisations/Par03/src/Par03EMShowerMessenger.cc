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
#include "Par03EMShowerModel.hh"
#include "Par03EMShowerMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

Par03EMShowerMessenger::Par03EMShowerMessenger(Par03EMShowerModel* aModel)
  : fModel(aModel)
{
  fDirectory = new G4UIdirectory("/Par03/fastSim/");
  fDirectory->SetGuidance(
    "Set mesh parameters for the example fast sim model.");

  fPrintCmd = new G4UIcmdWithoutParameter("/Par03/fastSim/print", this);
  fPrintCmd->SetGuidance("Print current settings.");

  fSigmaCmd = new G4UIcmdWithADoubleAndUnit(
    "/Par03/fastSim/transverseProfile/sigma", this);
  fSigmaCmd->SetGuidance("Set sigma parameter of 2D Gaussian distribution.");
  fSigmaCmd->SetParameterName("Sigma", false);
  fSigmaCmd->SetUnitCategory("Length");

  fAlphaCmd =
    new G4UIcmdWithADouble("/Par03/fastSim/longitudinalProfile/alpha", this);
  fAlphaCmd->SetGuidance("Set alpha parameter of Gamma distribution.");
  fAlphaCmd->SetParameterName("Alpha", false);

  fBetaCmd =
    new G4UIcmdWithADouble("/Par03/fastSim/longitudinalProfile/beta", this);
  fBetaCmd->SetGuidance("Set beta parameter of Gamma distribution.");
  fBetaCmd->SetParameterName("Beta", false);

  fNbOfHitsCmd = new G4UIcmdWithAnInteger("/Par03/fastSim/numberOfHits", this);
  fNbOfHitsCmd->SetGuidance(
    "Set number of (same energy) energy deposits created in fast simulation. "
    "Those deposits will be scored in the detector according to the readout of "
    "the sensitive detector.");
  fNbOfHitsCmd->SetParameterName("Number", false);

  fLongMaxDepthCmd =
    new G4UIcmdWithADouble("/Par03/fastSim/longitudinalProfile/maxDepth", this);
  fLongMaxDepthCmd->SetGuidance(
    "Set maximum shower depth used in parametrisation.");
  fLongMaxDepthCmd->SetGuidance("Expressed in units of radiation length.");
  fLongMaxDepthCmd->SetParameterName("Depth", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03EMShowerMessenger::~Par03EMShowerMessenger()
{
  delete fPrintCmd;
  delete fSigmaCmd;
  delete fAlphaCmd;
  delete fBetaCmd;
  delete fNbOfHitsCmd;
  delete fLongMaxDepthCmd;
  delete fDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03EMShowerMessenger::SetNewValue(G4UIcommand* aCommand,
                                         G4String aNewValues)
{
  if(aCommand == fPrintCmd)
  {
    fModel->Print();
  }
  else if(aCommand == fSigmaCmd)
  {
    fModel->SetSigma(fSigmaCmd->GetNewDoubleValue(aNewValues));
  }
  else if(aCommand == fAlphaCmd)
  {
    fModel->SetAlpha(fAlphaCmd->GetNewDoubleValue(aNewValues));
  }
  else if(aCommand == fBetaCmd)
  {
    fModel->SetBeta(fBetaCmd->GetNewDoubleValue(aNewValues));
  }
  else if(aCommand == fNbOfHitsCmd)
  {
    fModel->SetNbOfHits(fNbOfHitsCmd->GetNewIntValue(aNewValues));
  }
  else if(aCommand == fLongMaxDepthCmd)
  {
    fModel->SetLongMaxDepth(fLongMaxDepthCmd->GetNewDoubleValue(aNewValues));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par03EMShowerMessenger::GetCurrentValue(G4UIcommand* aCommand)
{
  G4String cv;

  if(aCommand == fSigmaCmd)
  {
    cv = fSigmaCmd->ConvertToString(fModel->GetSigma());
  }
  else if(aCommand == fAlphaCmd)
  {
    cv = fAlphaCmd->ConvertToString(fModel->GetAlpha());
  }
  else if(aCommand == fBetaCmd)
  {
    cv = fBetaCmd->ConvertToString(fModel->GetBeta());
  }
  else if(aCommand == fNbOfHitsCmd)
  {
    cv = fNbOfHitsCmd->ConvertToString(fModel->GetNbOfHits());
  }
  else if(aCommand == fLongMaxDepthCmd)
  {
    cv = fLongMaxDepthCmd->ConvertToString(fModel->GetLongMaxDepth());
  }
  return cv;
}
