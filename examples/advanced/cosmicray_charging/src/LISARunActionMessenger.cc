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

#include <sstream>

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
    std::istringstream is(t);
    is >> vl;
    runAction->SetAutoSeed(vl!=0);
  }

}

