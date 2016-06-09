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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
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
  colourNeutronCmd->AvailableForStates(G4State_Idle);

  colourGammaCmd = new G4UIcmdWithAString("/dmx/draw/gammaColour",this);
  colourGammaCmd->SetGuidance("Colour of gamma in the event");
  colourGammaCmd->SetGuidance("  Choice : white, grey, lgrey, black, red, green, blue, cyan(default), magenta, yellow, lgreen, lblue");
  colourGammaCmd->SetParameterName("colourGammaFlag", false);
  colourGammaCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourGammaCmd->AvailableForStates(G4State_Idle);

  colourOpticalCmd = new G4UIcmdWithAString("/dmx/draw/opticalColour",this);
  colourOpticalCmd->SetGuidance("Colour of gamma in the event");
  colourOpticalCmd->SetGuidance("  Choice : white(default), grey, lgrey, black, red, green, blue, cyan, magenta, yellow, lgreen, lblue");
  colourOpticalCmd->SetParameterName("colourOpticalFlag", false);
  colourOpticalCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourOpticalCmd->AvailableForStates(G4State_Idle);

  colourChargedPlusCmd = new G4UIcmdWithAString("/dmx/draw/chargedplusColour",this);
  colourChargedPlusCmd->SetGuidance("colour of chargedplus in the event");
  colourChargedPlusCmd->SetGuidance("  Choice : white, grey, lgrey, black, red (default), green, blue, cyan, magenta, yellow, lgreen, lblue");
  colourChargedPlusCmd->SetParameterName("ColourChargedPlusFlag", false);
  colourChargedPlusCmd->SetCandidates("white grey lgrey black red(default) green blue cyan magenta yellow lgreen lblue");
  colourChargedPlusCmd->AvailableForStates(G4State_Idle);

  colourChargedMinusCmd = new G4UIcmdWithAString("/dmx/draw/chargedminusColour",this);
  colourChargedMinusCmd->SetGuidance("colour of chargedminus in the event");
  colourChargedMinusCmd->SetGuidance("  Choice : white, grey, lgrey, black, red, green, blue(default), cyan, magenta, yellow, lgreen, lblue");
  colourChargedMinusCmd->SetParameterName("colourChargedMinusFlag", false);
  colourChargedMinusCmd->SetCandidates("white grey lgrey black red green blue cyan magenta yellow lgreen lblue");
  colourChargedMinusCmd->AvailableForStates(G4State_Idle);

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



