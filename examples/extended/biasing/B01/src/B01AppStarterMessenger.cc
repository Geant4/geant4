#include "B01AppStarterMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "B01VSimulation.hh"
#include "B01AppStarter.hh"

B01AppStarterMessenger::
B01AppStarterMessenger(B01AppStarter *appstarter)
  :
  fAppStarter(appstarter),
  fWeightRoulette(false),
  fSim(0),
  fApp(false),
  fVisAppComand(new G4UIcommand("/B01/app/vis", this)),
  fTimedAppComand(new G4UIcmdWithAnInteger("/B01/app/timed", this))
{
  if (!fVisAppComand) {
    G4std::G4Exception("B01AppStarterMessenger new failed!");
  }    
  if (!fTimedAppComand) {
    G4std::G4Exception("B01AppStarterMessenger new failed!");
  }
  fTimedAppComand->SetDefaultValue(0);
  
  G4String base("/B01/sim/");
  B01SimNameVec names = fSimfac.GetSimulationNames();
  for (B01SimNameVec::iterator it = names.begin();
       it != names.end(); ++it){
    G4String cmdname(base + *it);
    G4UIcmdWithAnInteger *cmd = new G4UIcmdWithAnInteger(cmdname, this);
    if (!cmd) {
      G4std::G4Exception("B01AppStarterMessenger new failed!");
    }
    fSimComands[cmd] = *it;
    cmd->SetDefaultValue(0);
  }
  
}

B01AppStarterMessenger::
~B01AppStarterMessenger(){
  delete fVisAppComand;
  delete fTimedAppComand;
  for (B01SimComands::iterator it = fSimComands.begin();
       it != fSimComands.end(); it++){
    delete it->first;
    fSimComands.erase(it);
  }
  
}
void B01AppStarterMessenger::SetNewValue(G4UIcommand* pCmd,
					 G4String szValue) {
  B01SimComands::iterator simit = fSimComands.end();
  for (B01SimComands::iterator it = fSimComands.begin();
       it != fSimComands.end(); ++it){
    if (it->first == pCmd) {
      simit=it;
      break;
    }
  }

  if( simit!=fSimComands.end()) {
    if (fSim) {
      G4std::G4cout << "B01AppStarterMessenger::SetNewValue: an simulation already exists, will be deleted" << G4endl;      
      delete fSim;
    }
    else {
      G4int i = simit->first->GetNewIntValue(szValue);
      if (i>0) {
	fWeightRoulette = true;
      }
      fSim = fSimfac.Create(simit->second);
      G4std::G4cout << "B01AppStarterMessenger::SetNewValue: creating: " 
	     << simit->second << G4endl;
      fSim->SetWeightRoulette(fWeightRoulette);
      fAppStarter->SetSimulation(fSim);
    }
  }
  
  if (pCmd==fVisAppComand){
    if (fApp) {
      G4std::G4cout << "B01AppStarterMessenger::SetNewValue: an application already exists" << G4endl;
    }
    else {
      fApp = true;
      fAppStarter->CreateVisApplication();
    }
  }
  if (pCmd==fTimedAppComand){
    if (fApp) {
      G4std::G4cout << "B01AppStarterMessenger::SetNewValue: an application already exists" << G4endl;
    }
    else {
      fApp = true;
      G4int time = fTimedAppComand->GetNewIntValue(szValue);
      fAppStarter->CreateTimedApplication(time);
    }
  }
  if (fApp && fSim) {
    fAppStarter->Run();
  }
  
}



