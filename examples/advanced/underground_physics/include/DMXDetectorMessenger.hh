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

#ifndef DMXDetectorMessenger_h
#define DMXDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DMXDetectorConstruction;
class G4UIdirectory;

class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class DMXDetectorMessenger: public G4UImessenger {

  public:
    DMXDetectorMessenger(DMXDetectorConstruction*);
   ~DMXDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

  private:
    DMXDetectorConstruction*  detectorConstruction;   
  
    G4UIcmdWithADoubleAndUnit* RoomEKineCutCmd;
    G4UIcmdWithADoubleAndUnit* EKineCutCmd;
    G4UIcmdWithADoubleAndUnit* RoomTimeCutCmd;
    G4UIcmdWithADoubleAndUnit* TimeCutCmd;
  //    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

