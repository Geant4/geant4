#include "Tst14RunMessenger.hh"

#include "Tst14RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

Tst14RunMessenger::Tst14RunMessenger(Tst14RunAction* RA)
:runAction (RA)
{
  plotDir = new G4UIdirectory("/plots/");
  plotDir->SetGuidance("plot control");

  sethistNameCmd = new G4UIcmdWithAString("/plots/sethistName",this);
  sethistNameCmd->SetGuidance("set name for the histogram file"); 

  setnbinStepCmd = new G4UIcmdWithAnInteger("/plots/setnbinStep",this);
  setnbinStepCmd->SetGuidance("set nb of bins in #step plot");
  setnbinStepCmd->SetParameterName("nbinStep",true);

  setSteplowCmd = new G4UIcmdWithADouble("/plots/setSteplow",this);
  setSteplowCmd->SetGuidance("set lower limit for #step plot ");
  setSteplowCmd->SetParameterName("Steplow",true);

  setStephighCmd = new G4UIcmdWithADouble("/plots/setStephigh",this);
  setStephighCmd->SetGuidance("set upper limit for #step plot ");
  setStephighCmd->SetParameterName("Stephigh",true);

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

  setnbinTbCmd = new G4UIcmdWithAnInteger("/plots/setnbinTb",this);
  setnbinTbCmd->SetGuidance("set nb of bins in Ebackscattering plot");
  setnbinTbCmd->SetParameterName("nbinTb",true);

  setTblowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTblow",this);
  setTblowCmd->SetGuidance("set lower limit for Ebackscattered plot ");
  setTblowCmd->SetParameterName("Tblow",true);
  setTblowCmd->SetDefaultUnit("MeV");

  setTbhighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setTbhigh",this);
  setTbhighCmd->SetGuidance("set upper limit for Ebackscattered plot ");
  setTbhighCmd->SetParameterName("Tbhigh",true);
  setTbhighCmd->SetDefaultUnit("MeV");

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

  setnbinRCmd = new G4UIcmdWithAnInteger("/plots/setnbinR",this);
  setnbinRCmd->SetGuidance("set nb of bins in R plot");
  setnbinRCmd->SetParameterName("nbinR",true);

  setRlowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setRlow",this);
  setRlowCmd->SetGuidance("set lower limit for R plot ");
  setRlowCmd->SetParameterName("Rlow",true);
  setRlowCmd->SetDefaultUnit("mm");

  setRhighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setRhigh",this);
  setRhighCmd->SetGuidance("set upper limit for R plot ");
  setRhighCmd->SetParameterName("Rhigh",true);
  setRhighCmd->SetDefaultUnit("mm");

  setnbinzvertexCmd = new G4UIcmdWithAnInteger("/plots/setnbinzvertex",this);
  setnbinzvertexCmd->SetGuidance("set nb of bins in Z vertex plot");
  setnbinzvertexCmd->SetParameterName("nbinZ",true);

  setzlowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setzlow",this);
  setzlowCmd->SetGuidance("set lower limit for Z vertex plot ");
  setzlowCmd->SetParameterName("zlow",true);
  setzlowCmd->SetDefaultUnit("mm");

  setzhighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setzhigh",this);
  setzhighCmd->SetGuidance("set upper limit for Z vertex plot ");
  setzhighCmd->SetParameterName("zhigh",true);
  setzhighCmd->SetDefaultUnit("mm");

  setnbinThCmd = new G4UIcmdWithAnInteger("/plots/setnbinTh",this);
  setnbinThCmd->SetGuidance("set nb of bins in Theta transmitted plot");
  setnbinThCmd->SetParameterName("nbinTh",true);

  setThlowCmd = new G4UIcmdWithADoubleAndUnit("/plots/setThlow",this);
  setThlowCmd->SetGuidance("set lower limit for Theta transmitted plot ");
  setThlowCmd->SetParameterName("Thlow",true);
  setThlowCmd->SetDefaultUnit("degree");

  setThhighCmd = new G4UIcmdWithADoubleAndUnit("/plots/setThhigh",this);
  setThhighCmd->SetGuidance("set upper limit for Theta transmitted plot ");
  setThhighCmd->SetParameterName("Thhigh",true);
  setThhighCmd->SetDefaultUnit("degree");

  setnbinThbackCmd = new G4UIcmdWithAnInteger("/plots/setnbinThback",this);
  setnbinThbackCmd->SetGuidance("set nb of bins in backscattering Theta plot");
  setnbinThbackCmd->SetParameterName("nbinThback",true);

  setThlowbackCmd = new G4UIcmdWithADoubleAndUnit("/plots/setThlowback",this);
  setThlowbackCmd->SetGuidance("set lower limit for backscattering Theta plot ");
  setThlowbackCmd->SetParameterName("Thlowback",true);
  setThlowbackCmd->SetDefaultUnit("degree");

  setThhighbackCmd = new G4UIcmdWithADoubleAndUnit("/plots/setThhighback",this);
  setThhighbackCmd->SetGuidance("set upper limit for backscattering Theta plot ");
  setThhighbackCmd->SetParameterName("Thhighback",true);
  setThhighbackCmd->SetDefaultUnit("degree");


}

Tst14RunMessenger::~Tst14RunMessenger()
{
  delete sethistNameCmd;

  delete setnbinStepCmd;
  delete setSteplowCmd;
  delete setStephighCmd;

  delete setnbinEnCmd;
  delete setEnlowCmd;
  delete setEnhighCmd;

  delete setnbinGammaCmd;
  delete setElowGammaCmd;
  delete setEhighGammaCmd;

  delete setnbinTtCmd;
  delete setTtlowCmd;
  delete setTthighCmd;

  delete setnbinTbCmd;
  delete setTblowCmd;
  delete setTbhighCmd;

  delete setnbinTsecCmd;
  delete setTseclowCmd;
  delete setTsechighCmd;

  delete setnbinRCmd;
  delete setRlowCmd;
  delete setRhighCmd;

  delete setnbinzvertexCmd;
  delete setzlowCmd;
  delete setzhighCmd;

  delete setnbinThCmd;
  delete setThlowCmd;
  delete setThhighCmd;

  delete setnbinThbackCmd;
  delete setThlowbackCmd;
  delete setThhighbackCmd;

  delete plotDir;

}

void Tst14RunMessenger::
         SetNewValue(G4UIcommand* command,G4String newValues)
{
  if( command == sethistNameCmd)
    runAction
    ->SethistName(newValues) ;
    
  if( command == setnbinStepCmd)
    runAction
    ->SetnbinStep(setnbinStepCmd->GetNewIntValue(newValues));

  if( command == setSteplowCmd)
    runAction
    ->SetSteplow( setSteplowCmd->GetNewDoubleValue(newValues));

  if( command == setStephighCmd)
    runAction
    ->SetStephigh( setStephighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinEnCmd)
    runAction
    ->SetnbinEn(setnbinEnCmd->GetNewIntValue(newValues));

  if( command == setEnlowCmd)
    runAction
    ->SetEnlow( setEnlowCmd->GetNewDoubleValue(newValues));

  if( command == setEnhighCmd)
    runAction
    ->SetEnhigh( setEnhighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinGammaCmd)
    runAction
    ->SetnbinGamma(setnbinGammaCmd->GetNewIntValue(newValues));

  if( command == setElowGammaCmd)
    runAction
    ->SetElowGamma( setElowGammaCmd->GetNewDoubleValue(newValues));

  if( command == setEhighGammaCmd)
    runAction
    ->SetEhighGamma( setEhighGammaCmd->GetNewDoubleValue(newValues));

  if( command == setnbinTtCmd)
    runAction
    ->SetnbinTt(setnbinTtCmd->GetNewIntValue(newValues));

  if( command == setTtlowCmd)
    runAction
    ->SetTtlow( setTtlowCmd->GetNewDoubleValue(newValues));

  if( command == setTthighCmd)
    runAction
    ->SetTthigh( setTthighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinTbCmd)
    runAction
    ->SetnbinTb(setnbinTbCmd->GetNewIntValue(newValues));

  if( command == setTblowCmd)
    runAction
    ->SetTblow( setTblowCmd->GetNewDoubleValue(newValues));

  if( command == setTbhighCmd)
    runAction
    ->SetTbhigh( setTbhighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinTsecCmd)
    runAction
    ->SetnbinTsec(setnbinTsecCmd->GetNewIntValue(newValues));

  if( command == setTseclowCmd)
    runAction
    ->SetTseclow( setTseclowCmd->GetNewDoubleValue(newValues));

  if( command == setTsechighCmd)
    runAction
    ->SetTsechigh( setTsechighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinRCmd)
    runAction
    ->SetnbinR(setnbinRCmd->GetNewIntValue(newValues));

  if( command == setRlowCmd)
    runAction
    ->SetRlow( setRlowCmd->GetNewDoubleValue(newValues));

  if( command == setRhighCmd)
    runAction
    ->SetRhigh( setRhighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinzvertexCmd)
    runAction
    ->Setnbinzvertex(setnbinzvertexCmd->GetNewIntValue(newValues));

  if( command == setzlowCmd)
    runAction
    ->Setzlow( setzlowCmd->GetNewDoubleValue(newValues));

  if( command == setzhighCmd)
    runAction
    ->Setzhigh( setzhighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinThCmd)
    runAction
    ->SetnbinTh(setnbinThCmd->GetNewIntValue(newValues));

  if( command == setThlowCmd)
    runAction
    ->SetThlow( setThlowCmd->GetNewDoubleValue(newValues));

  if( command == setThhighCmd)
    runAction
    ->SetThhigh( setThhighCmd->GetNewDoubleValue(newValues));

  if( command == setnbinThbackCmd)
    runAction
    ->SetnbinThBack(setnbinThbackCmd->GetNewIntValue(newValues));

  if( command == setThlowbackCmd)
    runAction
    ->SetThlowBack( setThlowbackCmd->GetNewDoubleValue(newValues));

  if( command == setThhighbackCmd)
    runAction
    ->SetThhighBack( setThhighbackCmd->GetNewDoubleValue(newValues));


}

   
