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
// RunActionMessenger header
// --------------------------------------------------------------

#ifndef DMXRunActionMessenger_h
#define DMXRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DMXRunAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXRunActionMessenger: public G4UImessenger
{
  public:
    DMXRunActionMessenger(DMXRunAction*);
   ~DMXRunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DMXRunAction*   DMXRun;

    G4UIcmdWithAString*   SaveHitsCmd;
    G4UIcmdWithAString*   SavePmtCmd;
    G4UIcmdWithAString*   SaveHistFileCmd;
};

#endif

