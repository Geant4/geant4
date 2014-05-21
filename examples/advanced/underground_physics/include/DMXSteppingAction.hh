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
    DMXSteppingAction();
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

    const DMXEventAction*  evtAction;  //pointer to event action


};

#endif
