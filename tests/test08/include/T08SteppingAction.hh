// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08SteppingAction.hh,v 1.1 1999-01-08 16:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08SteppingAction_h
#define T08SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class T08DetectorConstruction;
class T08EventAction;

class T08SteppingAction : public G4UserSteppingAction
{
  public:
    T08SteppingAction(T08DetectorConstruction* myDC,T08EventAction* myEA);
    ~T08SteppingAction(){};

    void UserSteppingAction();
    
  private:
    T08DetectorConstruction* myDetector;
    T08EventAction* eventAction;
};

#endif
