// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossMessenger.cc,v 1.3 2000-11-09 15:52:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossMessenger.hh"

#include "G4VEnergyLoss.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessenger::G4EnergyLossMessenger()
{
  eLossDirectory = new G4UIdirectory("/process/eLoss/");
  eLossDirectory->SetGuidance("Commands for G4VEnergyLoss.");
         
  RndmStepCmd = new G4UIcmdWithABool("/process/eLoss/rndmStep",this);
  RndmStepCmd->SetGuidance("Randomize the proposed step by eLoss.");
  RndmStepCmd->SetParameterName("choice",true);
  RndmStepCmd->SetDefaultValue(false);
  RndmStepCmd->AvailableForStates(Idle);
  
  EnlossFlucCmd = new G4UIcmdWithABool("/process/eLoss/fluct",this);
  EnlossFlucCmd->SetGuidance("Switch true/false the energy loss fluctuations.");
  EnlossFlucCmd->SetParameterName("choice",true);
  EnlossFlucCmd->SetDefaultValue(true);
  EnlossFlucCmd->AvailableForStates(Idle);

  SubSecCmd = new G4UIcmdWithABool("/process/eLoss/subsec",this);
  SubSecCmd->SetGuidance("Switch true/false the subcutoff generation.");
  SubSecCmd->SetParameterName("choice",true);
  SubSecCmd->SetDefaultValue(true);
  SubSecCmd->AvailableForStates(Idle);

  MinSubSecCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/minsubsec",this);
  MinSubSecCmd->SetGuidance("Set the min. cut for subcutoff delta in range.");
  MinSubSecCmd->SetParameterName("rcmin",true);
  MinSubSecCmd->AvailableForStates(Idle);

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
  StepFuncCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossMessenger::~G4EnergyLossMessenger()
{
  delete RndmStepCmd;
  delete EnlossFlucCmd;
  delete SubSecCmd;
  delete MinSubSecCmd;
  delete StepFuncCmd;
  delete eLossDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == RndmStepCmd)
   { G4VEnergyLoss::SetRndmStep(RndmStepCmd->GetNewBoolValue(newValue));
   }
   
  if (command == EnlossFlucCmd)
   { G4VEnergyLoss::SetEnlossFluc(EnlossFlucCmd->GetNewBoolValue(newValue));
   }

  if (command == SubSecCmd)
   { G4VEnergyLoss::SetSubSec(SubSecCmd->GetNewBoolValue(newValue));
   }

  if (command == MinSubSecCmd)
   { G4VEnergyLoss::SetMinDeltaCutInRange(MinSubSecCmd->GetNewDoubleValue(newValue));
   }

  if (command == StepFuncCmd)
   {
     G4double v1,v2;
     char unts[30];
     const char* t = newValue;
     G4std::istrstream is((char*)t);
     is >> v1 >> v2 >> unts;
     G4String unt = unts;
     v2 *= G4UIcommand::ValueOf(unt);
     G4VEnergyLoss::SetStepFunction(v1,v2);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
