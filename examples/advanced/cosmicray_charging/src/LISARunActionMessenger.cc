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

