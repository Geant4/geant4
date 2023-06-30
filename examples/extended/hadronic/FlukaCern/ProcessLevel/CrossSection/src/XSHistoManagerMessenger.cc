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
///  \file XSHistoManagerMessenger.cc
///  \brief UI commands for XS study.
//
//  Author: G.Hugo, 06 January 2023
//
// ***************************************************************************
//
//      XSHistoManagerMessenger
//
///  UI commands for XS study.
//
// ***************************************************************************

#include "XSHistoManagerMessenger.hh"

#include "XSHistoManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XSHistoManagerMessenger::XSHistoManagerMessenger(XSHistoManager* const histoManager)
  : fHisto(histoManager),
    fOutputFileNameCmd(G4UIcmdWithAString("/allXS/outputFileName", this)),
    fParticleNameCmd(G4UIcmdWithAString("/allXS/particleName", this)),
    fElementNameCmd(G4UIcmdWithAString("/allXS/elementName", this)),
    fNonElementaryMaterialNameCmd(G4UIcmdWithAString("/allXS/nonElementaryMaterialName", this)),
    fNumBinsCmd(G4UIcmdWithAnInteger("/allXS/numBins", this)),
    fMinKineticEnergyCmd(G4UIcmdWithADoubleAndUnit("/allXS/minKineticEnergy", this)),
    fMaxKineticEnergyCmd(G4UIcmdWithADoubleAndUnit("/allXS/maxKineticEnergy", this))
{
  fOutputFileNameCmd.SetGuidance("Set output file name (histograms).");
  
  fParticleNameCmd.SetGuidance("Set particle name.");
  fParticleNameCmd.SetParameterName("particleName", false);
  fParticleNameCmd.AvailableForStates(G4State_PreInit, G4State_Idle);
 
  fElementNameCmd.SetGuidance("Set target element name.");
  fElementNameCmd.SetParameterName("elementName", false);
  fElementNameCmd.AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fNonElementaryMaterialNameCmd.SetGuidance("Set target material name (in case not elementary).");
  fNonElementaryMaterialNameCmd.SetParameterName("nonElementaryMaterialName", false);
  fNonElementaryMaterialNameCmd.AvailableForStates(G4State_PreInit, G4State_Idle);

  fNumBinsCmd.SetGuidance("Set number of bins in kinetic energy.");
  fNumBinsCmd.SetParameterName("numBins", false);
  fNumBinsCmd.AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fMinKineticEnergyCmd.SetGuidance("Set min kinetic energy");
  fMinKineticEnergyCmd.SetParameterName("MinKineticEnergy", false);
  fMinKineticEnergyCmd.SetUnitCategory("Energy");
  fMinKineticEnergyCmd.AvailableForStates(G4State_PreInit, G4State_Idle);

  fMaxKineticEnergyCmd.SetGuidance("Set max kinetic energy");
  fMaxKineticEnergyCmd.SetParameterName("MaxKineticEnergy", false);
  fMaxKineticEnergyCmd.SetUnitCategory("Energy");
  fMaxKineticEnergyCmd.AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XSHistoManagerMessenger::SetNewValue(G4UIcommand* command, G4String value) {

  if (command == &fOutputFileNameCmd) {
    fHisto->SetOutputFileName(value);
  }
  else if (command == &fParticleNameCmd) {
    fHisto->SetParticle(value);
  }
  else if (command == &fElementNameCmd) {
    fHisto->SetElement(value);
  }
  else if (command == &fNonElementaryMaterialNameCmd) {
    fHisto->SetMaterial(value);
  }
  else if (command == &fNumBinsCmd) {
    fHisto->SetNumberOfBins(fNumBinsCmd.GetNewIntValue(value));
  }
  else if (command == &fMinKineticEnergyCmd) { 
    fHisto->SetMinKinEnergy(fMinKineticEnergyCmd.GetNewDoubleValue(value));
  }
  else if (command == &fMaxKineticEnergyCmd) { 
    fHisto->SetMaxKinEnergy(fMaxKineticEnergyCmd.GetNewDoubleValue(value));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
