#include "Tst33AppStarterMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "Tst33AppStarter.hh"

Tst33AppStarterMessenger::
Tst33AppStarterMessenger(Tst33AppStarter &appstarter)
  :
  G4UImessenger(),
  fAppStarter(appstarter),
  fMassGeoCmd(new G4UIcommand("/Tst33/MassGeometry", this)),
  fParallelGeoCmd(new G4UIcommand("/Tst33/ParallelGeometry",this)),
  fScoringCmd(new G4UIcommand("/Tst33/Scoring",this)),
  fImpCmd(new G4UIcommand("/Tst33/ImportanceSampling",this)),
  fWWRCmd(new G4UIcommand("/Tst33/WeightRoulette",this)),
  fClearSmaplingCmd(new G4UIcommand("/Tst33/ClearSampling",this)),
  fConfigureSamplingCmd(new G4UIcommand("/Tst33/ConfigureSampling",this)),
  fVisAppComand(new G4UIcommand("/Tst33/Visualization",this)),
  fTimedAppComand(new G4UIcmdWithAnInteger("/Tst33/Timed",this)),
  fPostRunCmd(new G4UIcommand("/Tst33/PostRun",this)),
  fRunCmd(new G4UIcmdWithAnInteger("/Tst33/Run", this))
{
}

Tst33AppStarterMessenger::
~Tst33AppStarterMessenger(){
  delete fMassGeoCmd;
  delete fParallelGeoCmd;
  delete fScoringCmd;
  delete fImpCmd;
  delete fWWRCmd;
  delete fClearSmaplingCmd;
  delete fConfigureSamplingCmd;
  delete fVisAppComand;
  delete fTimedAppComand;
  delete fPostRunCmd;
  delete fRunCmd;
}


void Tst33AppStarterMessenger::SetNewValue(G4UIcommand* pCmd,
					 G4String szValue) {
  if (pCmd==fMassGeoCmd) {
    fAppStarter.CreateMassGeometry();
  }
  if (pCmd==fParallelGeoCmd) {
    fAppStarter.CreateParallelGeometry();
  }
  if (pCmd==fScoringCmd) {
    fAppStarter.CreateScorer();
  }
  if (pCmd==fImpCmd) {
    fAppStarter.CreateIStore();
  }
  if (pCmd==fWWRCmd) {
    fAppStarter.CreateWeightRoulette();
  }
  if (pCmd==fClearSmaplingCmd) {
    fAppStarter.ClearSampling();
  }
  if (pCmd==fConfigureSamplingCmd) {
    fAppStarter.ConfigureSampling();
  }
  if (pCmd==fVisAppComand) {
    fAppStarter.CreateVisApplication();
  }
  if (pCmd==fTimedAppComand) {
    G4int time = fTimedAppComand->GetNewIntValue(szValue);
    fAppStarter.CreateTimedApplication(time);
  }
  if (pCmd==fPostRunCmd) {
    fAppStarter.PostRun();
  }
  if (pCmd==fRunCmd) {
    G4int nevents = fRunCmd->GetNewIntValue(szValue);
    fAppStarter.Run(nevents);
  }
  
}



