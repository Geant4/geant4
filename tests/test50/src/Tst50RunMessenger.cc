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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 

#include "Tst50RunMessenger.hh"
#include "Tst50RunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

Tst50RunMessenger::Tst50RunMessenger(Tst50RunAction* tst50run)
  :run(tst50run) 
{
  runDir = new  G4UIdirectory("/run/test/");
  runDir -> SetGuidance("run GUI");

  transmissionTestCmd = new G4UIcmdWithAString("/run/test/trans",this);
  transmissionTestCmd -> SetGuidance("if you want Transmission/Backscattering test on");
  transmissionTestCmd -> SetGuidance("if you want CSDARange/StoppingPower test off");
  transmissionTestCmd -> SetGuidance("Choice : on, off (default)");
  transmissionTestCmd -> SetParameterName("choice",true);
  transmissionTestCmd -> SetDefaultValue("off");
  transmissionTestCmd -> SetCandidates("on off");
  transmissionTestCmd -> AvailableForStates(G4State_Idle);
}

Tst50RunMessenger::~Tst50RunMessenger()
{
  delete  runDir;
  delete  transmissionTestCmd;
}
 
void Tst50RunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == transmissionTestCmd){run->SetTransmissionTest(newValue);}
}
