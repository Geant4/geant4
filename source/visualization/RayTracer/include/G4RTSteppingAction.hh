// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTSteppingAction.hh,v 1.2 2000-03-09 15:36:34 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//////////////////////
//G4RTSteppingAction
/////////////////////


#ifndef G4RTSteppingAction_h
#define G4RTSteppingAction_h 1


#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4RTSteppingAction : public G4UserSteppingAction
{
  public:
    G4RTSteppingAction();
    virtual ~G4RTSteppingAction(){;}

    virtual void UserSteppingAction(const G4Step*);

  private:
    G4bool ignoreTransparency;

  public:
    inline void SetIgnoreTransparency(G4bool val)
    { ignoreTransparency = val; }
    inline G4bool GetIgnoreTransparency() const
    { return ignoreTransparency; }
};

#endif
