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
//
//$Id: RemSimSteppingActionMessenger.cc,v 1.5 2006-06-29 16:24:27 gunter Exp $// GEANT4 tag $Name: not supported by cvs2svn $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//

#include "RemSimSteppingActionMessenger.hh"
#include "RemSimSteppingAction.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

RemSimSteppingActionMessenger::RemSimSteppingActionMessenger(RemSimSteppingAction* SA)
:steppingAction(SA)
{
  stepDirectory = new G4UIdirectory("/step/");
  stepDirectory -> SetGuidance("Step control command.");

  hadronicCmd = new G4UIcmdWithAString("/step/hadronicVerbose",this);
  hadronicCmd -> SetGuidance("Set the verbose level of hadronic processes");
  hadronicCmd -> SetParameterName("choice", true);
  hadronicCmd -> SetCandidates("On Off");
  hadronicCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 
}
RemSimSteppingActionMessenger::~RemSimSteppingActionMessenger()
{
  delete hadronicCmd;
  delete stepDirectory;
}

void RemSimSteppingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
 if(command == hadronicCmd) 
 {
   steppingAction -> SetHadronicAnalysis(newValue);
   G4cout<< " The stepping verbose is switched on" << G4endl; 
 }
}

