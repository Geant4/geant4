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
// SteppingActionMessenger header
// --------------------------------------------------------------

#ifndef DMXSteppingActionMessenger_h
#define DMXSteppingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DMXSteppingAction;
class G4UIdirectory;
class G4UIcmdWithAString;


class DMXSteppingActionMessenger: public G4UImessenger {

  public:
    DMXSteppingActionMessenger(DMXSteppingAction*);
   ~DMXSteppingActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

  private:
    DMXSteppingAction*  steppingAction;   
  
    G4UIcmdWithAString*   colourNeutronCmd;
    G4UIcmdWithAString*   colourGammaCmd;
    G4UIcmdWithAString*   colourOpticalCmd;
    G4UIcmdWithAString*   colourChargedPlusCmd;
    G4UIcmdWithAString*   colourChargedMinusCmd;

};

#endif

