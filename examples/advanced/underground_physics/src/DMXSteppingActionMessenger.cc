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
// SteppingActionMessenger program
// --------------------------------------------------------------

#include "DMXSteppingActionMessenger.hh"

#include "DMXSteppingAction.hh"
#include "DMXEventActionMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

DMXSteppingActionMessenger::DMXSteppingActionMessenger
   (DMXSteppingAction* SA):steppingAction(SA) {

  colourNeutronCmd = new G4UIcmdWithAString("/dmx/draw/neutronColour",this);
  colourNeutronCmd->SetGuidance("Colour of neutron in the event");
  colourNeutronCmd->SetGuidance("  Choice : white, grey, lgrey, black, red, green, blue, cyan, magenta(default), yellow, lgreen, lblue");
  colourNeutronCmd->SetParameterName("colourNeutronFlag", false);
  colourNeutronCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourNeutronCmd->AvailableForStates(Idle);

  colourGammaCmd = new G4UIcmdWithAString("/dmx/draw/gammaColour",this);
  colourGammaCmd->SetGuidance("Colour of gamma in the event");
  colourGammaCmd->SetGuidance("  Choice : white, grey, lgrey, black, red, green, blue, cyan(default), magenta, yellow, lgreen, lblue");
  colourGammaCmd->SetParameterName("colourGammaFlag", false);
  colourGammaCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourGammaCmd->AvailableForStates(Idle);

  colourOpticalCmd = new G4UIcmdWithAString("/dmx/draw/opticalColour",this);
  colourOpticalCmd->SetGuidance("Colour of gamma in the event");
  colourOpticalCmd->SetGuidance("  Choice : white(default), grey, lgrey, black, red, green, blue, cyan, magenta, yellow, lgreen, lblue");
  colourOpticalCmd->SetParameterName("colourOpticalFlag", false);
  colourOpticalCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourOpticalCmd->AvailableForStates(Idle);

  colourChargedPlusCmd = new G4UIcmdWithAString("/dmx/draw/chargedplusColour",this);
  colourChargedPlusCmd->SetGuidance("colour of chargedplus in the event");
  colourChargedPlusCmd->SetGuidance("  Choice : white, grey, lgrey, black, red (default), green, blue, cyan, magenta, yellow, lgreen, lblue");
  colourChargedPlusCmd->SetParameterName("ColourChargedPlusFlag", false);
  colourChargedPlusCmd->SetCandidates("white grey lgrey black red(default) green blue cyan magenta yellow lgreen lblue");
  colourChargedPlusCmd->AvailableForStates(Idle);

  colourChargedMinusCmd = new G4UIcmdWithAString("/dmx/draw/chargedminusColour",this);
  colourChargedMinusCmd->SetGuidance("colour of chargedminus in the event");
  colourChargedMinusCmd->SetGuidance("  Choice : white, grey, lgrey, black, red, green, blue(default), cyan, magenta, yellow, lgreen, lblue");
  colourChargedMinusCmd->SetParameterName("colourChargedMinusFlag", false);
  colourChargedMinusCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourChargedMinusCmd->AvailableForStates(Idle);

}


//ooooooooooooooooooooooooooooooooooooooooo
DMXSteppingActionMessenger::~DMXSteppingActionMessenger() {

  delete colourNeutronCmd;  
  delete colourGammaCmd;  
  delete colourOpticalCmd;  
  delete colourChargedPlusCmd;  
  delete colourChargedMinusCmd;  
}


//ooooooooooooooooooooooooooooooooooooooooo
void DMXSteppingActionMessenger::SetNewValue(G4UIcommand* command, 
  G4String newValue) {

  if(command == colourNeutronCmd)
    steppingAction->SetColourNeutronFlag(newValue);

  if(command == colourGammaCmd)
    steppingAction->SetColourGammaFlag(newValue);

  if(command == colourOpticalCmd)
    steppingAction->SetColourOpticalFlag(newValue);

  if(command == colourChargedPlusCmd)
    steppingAction->SetColourChargedPlusFlag(newValue);

  if(command == colourChargedMinusCmd)
    steppingAction->SetColourChargedMinusFlag(newValue);

}



