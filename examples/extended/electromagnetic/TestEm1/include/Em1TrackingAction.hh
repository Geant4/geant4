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
// $Id: Em1TrackingAction.hh,v 1.3 2001-07-11 09:57:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef Em1TrackingAction_h
#define Em1TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Em1RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1TrackingAction : public G4UserTrackingAction {

  public:  
    Em1TrackingAction(Em1RunAction*);
   ~Em1TrackingAction() {};
   
    void PostUserTrackingAction(const G4Track*);
    
  private:
    Em1RunAction* runAction;  
};

#endif
