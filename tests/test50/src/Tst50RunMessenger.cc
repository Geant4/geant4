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
