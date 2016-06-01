// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.hh,v 1.3 1998/10/09 14:31:06 japost Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 

#ifndef ExN02SteppingAction_h
#define ExN02SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class ExN02DetectorConstruction;
class ExN02EventAction;

class ExN02SteppingAction : public G4UserSteppingAction
{
  public:
    ExN02SteppingAction(ExN02DetectorConstruction* myDC,ExN02EventAction* myEA);
    ~ExN02SteppingAction(){};

    void UserSteppingAction();
    
  private:
    ExN02DetectorConstruction* myDetector;
    ExN02EventAction* eventAction;
};

#endif
