// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05SteppingAction.hh,v 1.1 1999-01-07 16:06:14 gunter Exp $
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
    ~ExN05SteppingAction(){};

    void UserSteppingAction();

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif
