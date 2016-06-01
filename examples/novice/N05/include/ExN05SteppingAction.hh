// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05SteppingAction.hh,v 2.1 1998/07/12 02:42:17 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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
