#include "Tst17RunMessenger.hh"

#include "Tst17RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

Tst17RunMessenger::Tst17RunMessenger(Tst17RunAction* RA)
:runAction (RA)
{
  plotDir = new G4UIdirectory("/plots/");
  plotDir->SetGuidance("plot control");

  sethistNameCmd = new G4UIcmdWithAString("/plots/sethistName",this);
  sethistNameCmd->SetGuidance("set name for the histogram file"); 

  setnbinEnCmd = new G4UIcmdWithAnInteger("/plots/setnbinEn",this);
  setnbinEnCmd->SetGuidance("set nb of bins in Edep plot");
  setnbinEnCmd->SetParameterName("nbinE",true);

  setEnlowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setEnlow",this);
  setEnlowCmd->SetGuidance("set lower limit for Edep plot ");
  setEnlowCmd->SetParameterName("Elow",true);
  setEnlowCmd->SetDefaultUnit("MeV");

  setEnhighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setEnhigh",this);
  setEnhighCmd->SetGuidance("set upper limit for Edep plot ");
  setEnhighCmd->SetParameterName("Ehigh",true);
  setEnhighCmd->SetDefaultUnit("MeV");

  setnbinGammaCmd = new G4UIcmdWithAnInteger("/plots/setnbinGamma",this);
  setnbinGammaCmd->SetGuidance("set nb of bins in gamma spectrum plot");
  setnbinGammaCmd->SetParameterName("nbinGamma",true);

  setElowGammaCmd = new G4UIcmdWithADoubleAndUnit("/plots/setElowGamma",this);
  setElowGammaCmd->SetGuidance("set lower limit for gamma spectrum plot ");
  setElowGammaCmd->SetParameterName("ElowGamma",true);
  setElowGammaCmd->SetDefaultUnit("MeV");

  setEhighGammaCmd = new G4UIcmdWithADoubleAndUnit("/plots/setEhighGamma",this);
  setEhighGammaCmd->SetGuidance("set upper limit for gamma spectrum plot ");
  setEhighGammaCmd->SetParameterName("EhighGamma",true);
  setEhighGammaCmd->SetDefaultUnit("MeV");

  setnbinTtCmd = new G4UIcmdWithAnInteger("/plots/setnbinTt",this);
  setnbinTtCmd->SetGuidance("set nb of bins in Etransmitted plot");
  setnbinTtCmd->SetParameterName("nbinTt",true);

  setTtlowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTtlow",this);
  setTtlowCmd->SetGuidance("set lower limit for Etransmitted plot ");
  setTtlowCmd->SetParameterName("Ttlow",true);
  setTtlowCmd->SetDefaultUnit("MeV");

  setTthighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTthigh",this);
  setTthighCmd->SetGuidance("set upper limit for Etransmitted plot ");
  setTthighCmd->SetParameterName("Tthigh",true);
  setTthighCmd->SetDefaultUnit("MeV");

  setnbinTsecCmd = new G4UIcmdWithAnInteger("/plots/setnbinTsec",this);
  setnbinTsecCmd->SetGuidance("set nb of bins in charged Tsecondary plot");
  setnbinTsecCmd->SetParameterName("nbinTsec",true);

  setTseclowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTseclow",this);
  setTseclowCmd->SetGuidance("set lower limit for charged Tsecondary plot ");
  setTseclowCmd->SetParameterName("Tseclow",true);
  setTseclowCmd->SetDefaultUnit("MeV");

  setTsechighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTsechigh",this);
  setTsechighCmd->SetGuidance("set upper limit for charged Tsecondary plot ");
  setTsechighCmd->SetParameterName("Tsechigh",true);
  setTsechighCmd->SetDefaultUnit("MeV");
}

Tst17RunMessenger::~Tst17RunMessenger()
{
  delete sethistNameCmd;

  delete setnbinEnCmd;
  delete setEnlowCmd;
  delete setEnhighCmd;

  delete setnbinGammaCmd;
  delete setElowGammaCmd;
  delete setEhighGammaCmd;

  delete setnbinTtCmd;
  delete setTtlowCmd;
  delete setTthighCmd;

  delete setnbinTsecCmd;
  delete setTseclowCmd;
  delete setTsechighCmd;

  delete plotDir;

}

void Tst17RunMessenger::SetNewValue(G4UIcommand* command,G4String newValues){

  if( command == sethistNameCmd){
    runAction->SethistName(newValues) ;
  }

  if( command == setnbinEnCmd){
    runAction->SetnbinEn(setnbinEnCmd->GetNewIntValue(newValues));
  }

  if( command == setEnlowCmd){
   runAction->SetEnlow( setEnlowCmd->GetNewDoubleValue(newValues));
  }

  if( command == setEnhighCmd){
    runAction->SetEnhigh( setEnhighCmd->GetNewDoubleValue(newValues));
  }

  if( command == setnbinGammaCmd){
    runAction->SetnbinGamma(setnbinGammaCmd->GetNewIntValue(newValues));
  }

  if( command == setElowGammaCmd){
   runAction->SetElowGamma( setElowGammaCmd->GetNewDoubleValue(newValues));
  }

  if( command == setEhighGammaCmd){
    runAction->SetEhighGamma( setEhighGammaCmd->GetNewDoubleValue(newValues));
  }

  if( command == setnbinTtCmd){
    runAction->SetnbinTt(setnbinTtCmd->GetNewIntValue(newValues));
  }

  if( command == setTtlowCmd){
    runAction->SetTtlow( setTtlowCmd->GetNewDoubleValue(newValues));
  }

  if( command == setTthighCmd){
    runAction->SetTthigh( setTthighCmd->GetNewDoubleValue(newValues));
  }

  if( command == setnbinTsecCmd){
    runAction->SetnbinTsec(setnbinTsecCmd->GetNewIntValue(newValues));
  }

  if( command == setTseclowCmd){
    runAction->SetTseclow( setTseclowCmd->GetNewDoubleValue(newValues));
  }

  if( command == setTsechighCmd){
    runAction->SetTsechigh( setTsechighCmd->GetNewDoubleValue(newValues));
  }
}

   









