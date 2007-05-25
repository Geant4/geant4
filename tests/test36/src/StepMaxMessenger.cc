#include "StepMaxMessenger.hh"

#include "StepMax.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxMessenger::StepMaxMessenger(StepMax* stepM)
:stepMax(stepM)
{
  StepMaxDir = new G4UIdirectory("/testem/stepMax/");
  StepMaxDir->SetGuidance("histograms control");
   
  StepMaxCmd = new G4UIcommand("/testem/stepMax/absorber",this);
  StepMaxCmd->SetGuidance("Set max allowed step length in absorber k");
  //
  G4UIparameter* k = new G4UIparameter("k",'i',false);
  k->SetGuidance("absorber number : from 1 to MaxHisto-1");
  k->SetParameterRange("k>0");
  StepMaxCmd->SetParameter(k);
  //    
  G4UIparameter* sMax = new G4UIparameter("sMax",'d',false);
  sMax->SetGuidance("stepMax, expressed in choosen unit");
  sMax->SetParameterRange("sMax>0.");
  StepMaxCmd->SetParameter(sMax);
  //    
  G4UIparameter* unit = new G4UIparameter("unit",'s',false);
  StepMaxCmd->SetParameter(unit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMaxMessenger::~StepMaxMessenger()
{
  delete StepMaxCmd;
  delete StepMaxDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{ 
  if (command == StepMaxCmd)
   { G4int k; G4double sMax; 
     G4String unts;
     std::istringstream is(newValues);
     is >> k >> sMax >> unts;
     G4String unit = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);  
     stepMax->SetStepMax(k,sMax*vUnit);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
