
#ifndef ExE02EventAction_h
#define ExE02EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ExE02EventAction : public G4UserEventAction
{
  public:
    ExE02EventAction();
    virtual ~ExE02EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
