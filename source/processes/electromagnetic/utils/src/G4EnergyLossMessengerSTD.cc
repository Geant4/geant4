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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4EnergyLossMessengerSTD
//
// Author:        
// 
// Creation date: 
//
// Modifications: 
//
// 03.01.2002 V.Ivanchenko update to new design
//
// Class Description: 
//

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossMessengerSTD.hh"

#include "G4LossTableManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessengerSTD::G4EnergyLossMessengerSTD()
{
  eLossDirectory = new G4UIdirectory("/process/eLoss/");
  eLossDirectory->SetGuidance("Commands for G4VEnergyLoss.");
         
  RndmStepCmd = new G4UIcmdWithABool("/process/eLoss/rndmStep",this);
  RndmStepCmd->SetGuidance("Randomize the proposed step by eLoss.");
  RndmStepCmd->SetParameterName("rndmStep",true);
  RndmStepCmd->SetDefaultValue(false);
  RndmStepCmd->AvailableForStates(PreInit,Idle);
  
  EnlossFlucCmd = new G4UIcmdWithABool("/process/eLoss/fluct",this);
  EnlossFlucCmd->SetGuidance("Switch true/false the energy loss fluctuations.");
  EnlossFlucCmd->SetParameterName("fluct",true);
  EnlossFlucCmd->SetDefaultValue(true);
  EnlossFlucCmd->AvailableForStates(PreInit,Idle);

  SubSecCmd = new G4UIcmdWithABool("/process/eLoss/subsec",this);
  SubSecCmd->SetGuidance("Switch true/false the subcutoff generation.");
  SubSecCmd->SetParameterName("subsec",true);
  SubSecCmd->SetDefaultValue(true);
  SubSecCmd->AvailableForStates(PreInit,Idle);

  IntegCmd = new G4UIcmdWithABool("/process/eLoss/integral",this);
  IntegCmd->SetGuidance("Switch true/false the integration of cross section over step.");
  IntegCmd->SetParameterName("integ",true);
  IntegCmd->SetDefaultValue(true);
  IntegCmd->AvailableForStates(PreInit,Idle);

  MinSubSecCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/minsubsec",this);
  MinSubSecCmd->SetGuidance("Set the min. cut for subcutoff delta in range.");
  MinSubSecCmd->SetParameterName("minsubsec",true);
  MinSubSecCmd->AvailableForStates(PreInit,Idle);

  MinEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/minKinEnergy",this);
  MinEnCmd->SetGuidance("Set the min kinetic energy");
  MinEnCmd->SetParameterName("emin",true);
  MinEnCmd->SetUnitCategory("Energy");
  MinEnCmd->AvailableForStates(PreInit,Idle);

  MaxEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/maxKinEnergy",this);
  MaxEnCmd->SetGuidance("Set the max kinetic energy");
  MaxEnCmd->SetParameterName("emax",true);
  MaxEnCmd->SetUnitCategory("Energy");
  MaxEnCmd->AvailableForStates(PreInit,Idle);

  StepFuncCmd = new G4UIcommand("/process/eLoss/StepFunction",this);
  StepFuncCmd->SetGuidance("Set the energy loss step limitation parameters.");
  StepFuncCmd->SetGuidance("  dRoverR   : max Range variation per step");
  StepFuncCmd->SetGuidance("  finalRange: range for final step");
  
  G4UIparameter* dRoverRPrm = new G4UIparameter("dRoverR",'d',false);
  dRoverRPrm->SetGuidance("max Range variation per step (fractional number)");
  dRoverRPrm->SetParameterRange("dRoverR>0. && dRoverR<=1.");
  StepFuncCmd->SetParameter(dRoverRPrm);
  
  G4UIparameter* finalRangePrm = new G4UIparameter("finalRange",'d',false);
  finalRangePrm->SetGuidance("range for final step");
  finalRangePrm->SetParameterRange("finalRange>0.");
  StepFuncCmd->SetParameter(finalRangePrm);
  
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',true);
  unitPrm->SetGuidance("unit of finalRange");
  unitPrm->SetDefaultValue("mm");
  G4String unitCandidates = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitCandidates);
  
  StepFuncCmd->SetParameter(unitPrm);
  StepFuncCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessengerSTD::~G4EnergyLossMessengerSTD()
{
  delete RndmStepCmd;
  delete EnlossFlucCmd;  
  delete SubSecCmd;
  delete MinSubSecCmd;
  delete MinEnCmd;
  delete MaxEnCmd;
  delete StepFuncCmd;
  delete IntegCmd;
  delete eLossDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossMessengerSTD::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  G4LossTableManager* lossTables = G4LossTableManager::Instance();

  if (command == RndmStepCmd) {
    lossTables->SetRandomStep(RndmStepCmd->GetNewBoolValue(newValue));
  }

  if (command == IntegCmd) {
    lossTables->SetIntegral(IntegCmd->GetNewBoolValue(newValue));
  }
   
  if (command == EnlossFlucCmd) {
    lossTables->SetLossFluctuations(EnlossFlucCmd->GetNewBoolValue(newValue));
  }

  if (command == SubSecCmd) {
    lossTables->SetSubCutoff(SubSecCmd->GetNewBoolValue(newValue));
  }

  if (command == MinEnCmd) {
    lossTables->SetMinEnergy(MinEnCmd->GetNewDoubleValue(newValue));
  }

  if (command == MaxEnCmd) {
    lossTables->SetMaxEnergy(MaxEnCmd->GetNewDoubleValue(newValue));
  }

  if (command == StepFuncCmd) {
     G4double v1,v2;
     char unts[30];
     const char* t = newValue;
     G4std::istrstream is((char*)t);
     is >> v1 >> v2 >> unts;
     G4String unt = unts;
     v2 *= G4UIcommand::ValueOf(unt);
     lossTables->SetStepLimits(v1,v2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
