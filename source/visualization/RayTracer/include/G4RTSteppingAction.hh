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
