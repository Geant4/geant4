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
// $Id: ConfigDataMessenger.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - messenger for ConfigData
//
#include "ConfigDataMessenger.hh"

#include "ConfigData.hh"
#include "Configuration.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

ConfigDataMessenger::ConfigDataMessenger()
{  
  listDir = new G4UIdirectory("/Brem/");
  
  fBe_15pt8MeV_1_10_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Be_15pt18MeV_1_10_degrees",this);  
  fBe_15pt8MeV_1_10_degrees->SetGuidance("Configuration for Be_15pt18MeV_1_10_degrees");
  fBe_15pt8MeV_1_10_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fBe_15pt18MeV_30_90_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Be_15pt18MeV_30_90_degrees",this);  
  fBe_15pt18MeV_30_90_degrees->SetGuidance("Configuration for Be_15pt18MeV_30_90_degrees");
  fBe_15pt18MeV_30_90_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_15pt8MeV_1_10_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_15pt18MeV_1_10_degrees",this);  
  fAl_15pt8MeV_1_10_degrees->SetGuidance("Configuration for Al_15pt18MeV_1_10_degrees");
  fAl_15pt8MeV_1_10_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_15pt18MeV_30_90_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_15pt18MeV_30_90_degrees",this);  
  fAl_15pt18MeV_30_90_degrees->SetGuidance("Configuration for Al_15pt18MeV_30_90_degrees");
  fAl_15pt18MeV_30_90_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_15pt8MeV_1_10_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_15pt18MeV_1_10_degrees",this);  
  fPb_15pt8MeV_1_10_degrees->SetGuidance("Configuration for Pb_15pt18MeV_1_10_degrees");
  fPb_15pt8MeV_1_10_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_15pt18MeV_30_90_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_15pt18MeV_30_90_degrees",this);  
  fPb_15pt18MeV_30_90_degrees->SetGuidance("Configuration for Pb_15pt18MeV_30_90_degrees");
  fPb_15pt18MeV_30_90_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ///////////////////////////////
  fPb_10pt09MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_10pt09MeV_0_degrees",this);  
  fPb_10pt09MeV_0_degrees->SetGuidance("Configuration for Pb_10pt09MeV_0_degrees");
  fPb_10pt09MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_15pt18MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_15pt18MeV_0_degrees",this);  
  fPb_15pt18MeV_0_degrees->SetGuidance("Configuration for Pb_15pt18MeV_0_degrees");
  fPb_15pt18MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_20pt28MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_20pt28MeV_0_degrees",this);  
  fPb_20pt28MeV_0_degrees->SetGuidance("Configuration for Pb_20pt28MeV_0_degrees");
  fPb_20pt28MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_25pt38MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_25pt38MeV_0_degrees",this);  
  fPb_25pt38MeV_0_degrees->SetGuidance("Configuration for Pb_25pt38MeV_0_degrees");
  fPb_25pt38MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPb_30pt45MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Pb_30pt45MeV_0_degrees",this);  
  fPb_30pt45MeV_0_degrees->SetGuidance("Configuration for Pb_30pt45MeV_0_degrees");
  fPb_30pt45MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ///////////////////////////////
  fAl_10pt09MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_10pt09MeV_0_degrees",this);  
  fAl_10pt09MeV_0_degrees->SetGuidance("Configuration for Al_10pt09MeV_0_degrees");
  fAl_10pt09MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_15pt18MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_15pt18MeV_0_degrees",this);  
  fAl_15pt18MeV_0_degrees->SetGuidance("Configuration for Al_15pt18MeV_0_degrees");
  fAl_15pt18MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_20pt28MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_20pt28MeV_0_degrees",this);  
  fAl_20pt28MeV_0_degrees->SetGuidance("Configuration for Al_20pt28MeV_0_degrees");
  fAl_20pt28MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_25pt38MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_25pt38MeV_0_degrees",this);  
  fAl_25pt38MeV_0_degrees->SetGuidance("Configuration for Al_25pt38MeV_0_degrees");
  fAl_25pt38MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAl_30pt45MeV_0_degrees = new G4UIcmdWithoutParameter("/Brem/Geometry/Al_30pt45MeV_0_degrees",this);  
  fAl_30pt45MeV_0_degrees->SetGuidance("Configuration for Al_30pt45MeV_0_degrees");
  fAl_30pt45MeV_0_degrees->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fVerbose = new G4UIcmdWithoutParameter("/Brem/Verbose",this);  
  fVerbose->SetGuidance("Verbose");
  fVerbose->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Min radius of scoring volume
  scorerMinRadiusCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMinRadius",this);  
  scorerMinRadiusCmd->SetGuidance("Min radius of scoring volume");
  scorerMinRadiusCmd->SetParameterName("scorerMinRadius",false);
  scorerMinRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Max radius of scoring volume
  scorerMaxRadiusCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMaxRadius",this);  
  scorerMaxRadiusCmd->SetGuidance("Max radius of scoring volume");
  scorerMaxRadiusCmd->SetParameterName("scorerMaxRadius",false);
  scorerMaxRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Min theta of scoring volume
  scorerMinThetaCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMinTheta",this);  
  scorerMinThetaCmd->SetGuidance("Min theta of scoring volume");
  scorerMinThetaCmd->SetParameterName("scorerMinTheta",false);
  scorerMinThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Max theta of scoring volume
  scorerMaxThetaCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMaxTheta",this);  
  scorerMaxThetaCmd->SetGuidance("Max theta of scoring volume");
  scorerMaxThetaCmd->SetParameterName("scorerMaxTheta",false);
  scorerMaxThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Delta theta of scoring volume
  scorerDeltaThetaCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerDeltaTheta",this);  
  scorerDeltaThetaCmd->SetGuidance("Delta theta of scoring volume");
  scorerDeltaThetaCmd->SetParameterName("scorerDeltaTheta",false);
  scorerDeltaThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Min energy of scorers
  scorerMinEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMinEnergy",this);  
  scorerMinEnergyCmd->SetGuidance("Min energy of scorers");
  scorerMinEnergyCmd->SetParameterName("scorerMinEnergy",false);
  scorerMinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Max energy of scorers
  scorerMaxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerMaxEnergy",this);  
  scorerMaxEnergyCmd->SetGuidance("Max energy of scorers");
  scorerMaxEnergyCmd->SetParameterName("scorerMaxEnergy",false);
  scorerMaxEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Delta energy of scorers
  scorerDeltaEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Brem/Scoring/ScorerDeltaEnergy",this);  
  scorerDeltaEnergyCmd->SetGuidance("Delta energy of scorers");
  scorerDeltaEnergyCmd->SetParameterName("scorerDeltaEnergy",false);
  scorerDeltaEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Brem splitting
  bremSplittingActivation = new G4UIcmdWithABool("/Brem/BremSplitting/Active",this);  
  bremSplittingActivation->SetGuidance("Activate or deactivate brem splitting");
  bremSplittingActivation->SetParameterName("active",false);
  bremSplittingActivation->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  bremSplittingFactor = new G4UIcmdWithAnInteger("/Brem/BremSplitting/Factor",this);  
  bremSplittingFactor->SetGuidance("Splitting parameter - how many photons produced");
  bremSplittingFactor->SetParameterName("nSplit",false);
  bremSplittingFactor->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ConfigDataMessenger::~ConfigDataMessenger()
{
  // Cleanup
  delete fPb_10pt09MeV_0_degrees;
  delete fPb_15pt18MeV_0_degrees;
  delete fPb_20pt28MeV_0_degrees;
  delete fPb_25pt38MeV_0_degrees;
  delete fPb_30pt45MeV_0_degrees;
  
  delete fAl_10pt09MeV_0_degrees;
  delete fAl_15pt18MeV_0_degrees;
  delete fAl_20pt28MeV_0_degrees;
  delete fAl_25pt38MeV_0_degrees;

  delete listDir;

  delete scorerMinRadiusCmd;
  delete scorerMaxRadiusCmd;
  delete scorerMinThetaCmd;
  delete scorerMaxThetaCmd;
  delete scorerDeltaThetaCmd;

  delete scorerMinEnergyCmd;
  delete scorerMaxEnergyCmd;
  delete scorerDeltaEnergyCmd;

  delete fVerbose;

  delete bremSplittingActivation;
  delete bremSplittingFactor;
}

void ConfigDataMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fBe_15pt8MeV_1_10_degrees) { 
    Be_15pt18MeV_1_10_degrees::Initialise();
  } else if (command == fBe_15pt18MeV_30_90_degrees) { 
    Be_15pt18MeV_30_90_degrees::Initialise();
  }
  else if (command == fAl_15pt8MeV_1_10_degrees) { 
    Al_15pt18MeV_1_10_degrees::Initialise();
  } else if (command == fAl_15pt18MeV_30_90_degrees) { 
   Al_15pt18MeV_30_90_degrees::Initialise();
  }
  else if (command == fPb_15pt8MeV_1_10_degrees) { 
    Pb_15pt18MeV_1_10_degrees::Initialise();
  } else if (command == fPb_15pt18MeV_30_90_degrees) { 
    Pb_15pt18MeV_30_90_degrees::Initialise();
  } 
  else if (command ==  fPb_10pt09MeV_0_degrees) {
    Pb_10pt09MeV_0_degrees::Initialise();
  }
  else if (command ==  fPb_15pt18MeV_0_degrees) {
    Pb_15pt18MeV_0_degrees::Initialise();
  }
  else if (command ==  fPb_20pt28MeV_0_degrees) {
    Pb_20pt28MeV_0_degrees::Initialise();
  }
  else if (command ==  fPb_25pt38MeV_0_degrees) {
    Pb_25pt38MeV_0_degrees::Initialise();
  }
  else if (command ==  fPb_30pt45MeV_0_degrees) {
    Pb_30pt45MeV_0_degrees::Initialise();
  }
  else if (command ==  fAl_10pt09MeV_0_degrees) {
    Al_10pt09MeV_0_degrees::Initialise();
  }
  else if (command ==  fAl_15pt18MeV_0_degrees) {
    Al_15pt18MeV_0_degrees::Initialise();
  }
  else if (command ==  fAl_20pt28MeV_0_degrees) {
    Al_20pt28MeV_0_degrees::Initialise();
  }
  else if (command ==  fAl_25pt38MeV_0_degrees) {
    Al_25pt38MeV_0_degrees::Initialise();
  }
  else if (command ==  fAl_30pt45MeV_0_degrees) {
    Al_30pt45MeV_0_degrees::Initialise();
  } else if (command == scorerMinRadiusCmd) { 
    ConfigData::SetScorerMinRadius(scorerMinRadiusCmd->GetNewDoubleValue(newValue)); 
  } else if (command == scorerMaxRadiusCmd) {
    ConfigData::SetScorerMaxRadius(scorerMaxRadiusCmd->GetNewDoubleValue(newValue));  
  } else if (command == scorerMinThetaCmd) {
    ConfigData::SetScorerMinTheta(scorerMinThetaCmd->GetNewDoubleValue(newValue));   
  } else if (command == scorerMaxThetaCmd) { 
    ConfigData::SetScorerMaxTheta(scorerMaxThetaCmd->GetNewDoubleValue(newValue));  
  } else if (command == scorerDeltaThetaCmd) { 
    ConfigData::SetScorerDeltaTheta(scorerDeltaThetaCmd->GetNewDoubleValue(newValue));   
  } else if (command == scorerMinEnergyCmd) {
    ConfigData::SetScorerMinEnergy(scorerMinEnergyCmd->GetNewDoubleValue(newValue));   
  } else if (command == scorerMaxEnergyCmd) { 
    ConfigData::SetScorerMaxEnergy(scorerMaxEnergyCmd->GetNewDoubleValue(newValue));  
  } else if (command == scorerDeltaEnergyCmd) { 
    ConfigData::SetScorerDeltaEnergy(scorerDeltaEnergyCmd->GetNewDoubleValue(newValue));         
  } else if (command == fVerbose) { 
    ConfigData::SetVerbose(true);       
  } else if (command == bremSplittingActivation) {
    ConfigData::BremSplitting::SetActivation(bremSplittingActivation->GetNewBoolValue(newValue));
  } else if (command == bremSplittingFactor) {
    ConfigData::BremSplitting::SetFactor(bremSplittingFactor->GetNewIntValue(newValue));  
  }else {
   G4cerr << "ConfigDataMessenger: Command not found. " << newValue << G4endl;
 }
}

