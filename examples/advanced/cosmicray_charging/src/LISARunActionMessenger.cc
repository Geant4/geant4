// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISARunActionMessenger class                                     *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#include "LISARunActionMessenger.hh"

#include "LISARunAction.hh"


LISARunActionMessenger::LISARunActionMessenger
  (LISARunAction* runAct) : runAction(runAct){

  SetAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  SetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  SetAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  SetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  SetAutoSeedCmd->SetGuidance("Default = false");
  SetAutoSeedCmd->SetParameterName("autoSeed", false);
  SetAutoSeedCmd->AvailableForStates(G4State_Idle);

}


LISARunActionMessenger::~LISARunActionMessenger() {

  delete SetAutoSeedCmd;  
}


void LISARunActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) { 

  if(command == SetAutoSeedCmd) {
    G4int vl;
    const char* t = newValue;
    std::istrstream is((char*)t);
    is >> vl;
    runAction->SetAutoSeed(vl!=0);
  }

}

