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
// EventActionMessenger header
// --------------------------------------------------------------

#ifndef DMXEventActionMessenger_h
#define DMXEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DMXEventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;


class DMXEventActionMessenger: public G4UImessenger {

  public:
    DMXEventActionMessenger(DMXEventAction*);
   ~DMXEventActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DMXEventAction*     eventAction;   
  
    G4UIdirectory*        dmxDirectory;
    G4UIdirectory*        drawDirectory;
    G4UIcmdWithAString*   DrawTrksCmd;
    G4UIcmdWithAString*   DrawColsCmd;
    G4UIcmdWithABool*     DrawHitsCmd;    
    G4UIcmdWithABool*     SavePmtCmd;    
    G4UIcmdWithABool*     SaveHitsCmd;    
    G4UIcmdWithAnInteger* PrintCmd;    

};

#endif
