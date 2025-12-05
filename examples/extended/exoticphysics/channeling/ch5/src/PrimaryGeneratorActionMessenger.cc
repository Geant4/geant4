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
// gpaterno, October 2025
//
/// \file PrimaryGeneratorActionMessenger.cc
/// \brief Implementation of the PrimaryGeneratorActionMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* pg):
fPrimaryGeneratorAction(pg)
{     
    fUseGPSCmd = new G4UIcmdWithABool("/gun/UseGPS", this);
    fUseGPSCmd->SetGuidance("Choice : (1) true, (0) false.");
    fUseGPSCmd->SetParameterName("setGPS", true);
    fUseGPSCmd->SetDefaultValue(false);
    
         
    fPrimaryTypeCmd = new G4UIcmdWithAString("/gun/primaryType",this);
    fPrimaryTypeCmd->SetParameterName("primaryType",true);
    fPrimaryTypeCmd->SetDefaultValue("e-");
      
    fPrimaryEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryEnergy",this);
    fPrimaryEnergyCmd->SetParameterName("primaryEnergy",true);
    fPrimaryEnergyCmd->SetDefaultValue(6);
    fPrimaryEnergyCmd->SetDefaultUnit("GeV");
      
    fPrimaryRelSigmaEnergyCmd = new G4UIcmdWithADouble("/gun/primaryRelSigmaEnergy",this);
    fPrimaryRelSigmaEnergyCmd->SetParameterName("primaryRelSigmaEnergy",true);
    fPrimaryRelSigmaEnergyCmd->SetDefaultValue(1.0e-3);
    
    fPrimaryXCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryX",this);
    fPrimaryXCmd->SetParameterName("primaryX",true);
    fPrimaryXCmd->SetDefaultValue(0);
    fPrimaryXCmd->SetDefaultUnit("mm");

    fPrimaryYCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryY",this);
    fPrimaryYCmd->SetParameterName("primaryY",true);
    fPrimaryYCmd->SetDefaultValue(0);
    fPrimaryYCmd->SetDefaultUnit("mm");
  
    fPrimaryZCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryZ",this);
    fPrimaryZCmd->SetParameterName("primaryZ",true);
    fPrimaryZCmd->SetDefaultValue(0);
    fPrimaryZCmd->SetDefaultUnit("mm");

    fPrimaryTCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryT",this);
    fPrimaryTCmd->SetParameterName("primaryT",true);
    fPrimaryTCmd->SetDefaultValue(0);
    fPrimaryTCmd->SetDefaultUnit("ns");
      
    fPrimaryXpCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryXp",this);
    fPrimaryXpCmd->SetParameterName("primaryXp",true);
    fPrimaryXpCmd->SetDefaultValue(0);
    fPrimaryXpCmd->SetDefaultUnit("rad");
    
    fPrimaryYpCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryYp",this);
    fPrimaryYpCmd->SetParameterName("primaryYp",true);
    fPrimaryYpCmd->SetDefaultValue(0);
    fPrimaryYpCmd->SetDefaultUnit("rad");   
           
    fPrimarySigmaXCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaX",this);
    fPrimarySigmaXCmd->SetParameterName("primarySigmaX",true);
    fPrimarySigmaXCmd->SetDefaultValue(0);
    fPrimarySigmaXCmd->SetDefaultUnit("mm");

    fPrimarySigmaYCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaY",this);
    fPrimarySigmaYCmd->SetParameterName("primarySigmaY",true);
    fPrimarySigmaYCmd->SetDefaultValue(0);
    fPrimarySigmaYCmd->SetDefaultUnit("mm");

    fPrimarySigmaZCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaZ",this);
    fPrimarySigmaZCmd->SetParameterName("primarySigmaZ",true);
    fPrimarySigmaZCmd->SetDefaultValue(0);
    fPrimarySigmaZCmd->SetDefaultUnit("mm");
  
    fPrimarySigmaTCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaT",this);
    fPrimarySigmaTCmd->SetParameterName("primarySigmaT",true);
    fPrimarySigmaTCmd->SetDefaultValue(0);
    fPrimarySigmaTCmd->SetDefaultUnit("ns");
  
    fPrimarySigmaXpCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaXp",this);
    fPrimarySigmaXpCmd->SetParameterName("primarySigmaXp",true);
    fPrimarySigmaXpCmd->SetDefaultValue(1.0e-5);
    fPrimarySigmaXpCmd->SetDefaultUnit("rad");
  
    fPrimarySigmaYpCmd = new G4UIcmdWithADoubleAndUnit("/gun/primarySigmaYp",this);
    fPrimarySigmaYpCmd->SetParameterName("primarySigmaYp",true);
    fPrimarySigmaYpCmd->SetDefaultValue(1.0e-5);
    fPrimarySigmaYpCmd->SetDefaultUnit("rad");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
    delete fUseGPSCmd;
      
    delete fPrimaryTypeCmd; 
    delete fPrimaryEnergyCmd; 
    delete fPrimaryRelSigmaEnergyCmd;
    delete fPrimaryXCmd;  
    delete fPrimaryYCmd;
    delete fPrimaryZCmd;
    delete fPrimaryTCmd;
    delete fPrimaryXpCmd;
    delete fPrimaryYpCmd;   
    delete fPrimarySigmaXCmd; 
    delete fPrimarySigmaYCmd; 
    delete fPrimarySigmaZCmd;
    delete fPrimarySigmaTCmd; 
    delete fPrimarySigmaXpCmd; 
    delete fPrimarySigmaYpCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{                       
    if (command == fPrimaryTypeCmd) {
        fPrimaryGeneratorAction->SetType(newValue);
    }          
    if (command == fPrimaryEnergyCmd) {
        fPrimaryGeneratorAction->SetEnergy(fPrimaryEnergyCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimaryRelSigmaEnergyCmd) {
        fPrimaryGeneratorAction->
            SetRelSigmaEnergy(fPrimaryRelSigmaEnergyCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimaryXCmd) {
        fPrimaryGeneratorAction->SetX(fPrimaryXCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimaryYCmd) {
        fPrimaryGeneratorAction->SetY(fPrimaryYCmd->GetNewDoubleValue(newValue));
    }  
    if (command == fPrimaryZCmd) {
        fPrimaryGeneratorAction->SetZ(fPrimaryZCmd->GetNewDoubleValue(newValue));
    }  
    if (command == fPrimaryTCmd) {
        fPrimaryGeneratorAction->SetT(fPrimaryTCmd->GetNewDoubleValue(newValue));
    }  
    if (command == fPrimaryXpCmd) {
        fPrimaryGeneratorAction->SetXp(fPrimaryXpCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimaryYpCmd) {
        fPrimaryGeneratorAction->SetYp(fPrimaryYpCmd->GetNewDoubleValue(newValue));
    }  
    if (command == fPrimarySigmaXCmd) {
        fPrimaryGeneratorAction->SetSigmaX(fPrimarySigmaXCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimarySigmaYCmd) {
        fPrimaryGeneratorAction->SetSigmaY(fPrimarySigmaYCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimarySigmaZCmd) {
        fPrimaryGeneratorAction->SetSigmaZ(fPrimarySigmaZCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimarySigmaTCmd) {
        fPrimaryGeneratorAction->SetSigmaT(fPrimarySigmaTCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimarySigmaXpCmd) {
        fPrimaryGeneratorAction->SetSigmaXp(fPrimarySigmaXpCmd->GetNewDoubleValue(newValue));
    }
    if (command == fPrimarySigmaYpCmd) {
        fPrimaryGeneratorAction->SetSigmaYp(fPrimarySigmaYpCmd->GetNewDoubleValue(newValue));
    }       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

