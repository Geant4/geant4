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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// StackingActionMessenger program
// --------------------------------------------------------------

#include "DMXStackingActionMessenger.hh"
#include "DMXStackingAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "globals.hh"


DMXStackingActionMessenger::DMXStackingActionMessenger
  (DMXStackingAction* msa) :DMXAction(msa) {

  KillGammasCmd = new G4UIcmdWithABool("/dmx/KillGammasInConcrete",this);
  KillGammasCmd->
    SetGuidance("Kills gammas produced by neutrons in the concrete wall");
  KillGammasCmd->SetGuidance("Default = false");
  KillGammasCmd->SetParameterName("KillGammasFlag", false);
  KillGammasCmd->SetDefaultValue(false);

}

DMXStackingActionMessenger::~DMXStackingActionMessenger() {

  delete KillGammasCmd;
}

void DMXStackingActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) {

  if(command==KillGammasCmd) {
    G4int vl;
    const char* t = newValue;
    G4std::istrstream is((char*)t);
    is >> vl;
    DMXAction->SetKillGammasFlag(vl!=0);
  }

}

