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
// SteppingAction header
// --------------------------------------------------------------

#ifndef DMXSteppingAction_h
#define DMXSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"


class DMXSteppingActionMessenger;
class DMXEventAction;
class DMXAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXSteppingAction : public G4UserSteppingAction
{
  public:
    DMXSteppingAction(DMXEventAction*);
    virtual ~DMXSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  public:
    void SetColourNeutronFlag(G4String val)      {colourNeutronFlag     =val;};
    void SetColourGammaFlag(G4String val)        {colourGammaFlag       =val;};
    void SetColourOpticalFlag(G4String val)      {colourOpticalFlag     =val;};
    void SetColourChargedPlusFlag(G4String val)  {colourChargedPlusFlag =val;};
    void SetColourChargedMinusFlag(G4String val) {colourChargedMinusFlag=val;};

  private:
    G4String                    colourNeutronFlag;
    G4String                    colourGammaFlag;
    G4String                    colourOpticalFlag;
    G4String                    colourChargedPlusFlag;
    G4String                    colourChargedMinusFlag;

    DMXSteppingActionMessenger*  steppingMessenger;

    DMXEventAction*    evtAction;  //pointer to event action


};

#endif
