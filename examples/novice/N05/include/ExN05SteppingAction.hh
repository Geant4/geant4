// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05SteppingAction.hh,v 1.2 1999-04-16 12:05:00 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05SteppingAction_h
#define ExN05SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class ExN05SteppingAction : public G4UserSteppingAction
{
  public:
    ExN05SteppingAction();
    virtual ~ExN05SteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif
