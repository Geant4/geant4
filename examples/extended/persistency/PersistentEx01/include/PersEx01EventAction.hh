
#ifndef PersEx01EventAction_h
#define PersEx01EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class PersEx01EventAction : public G4UserEventAction
{
  public:
    PersEx01EventAction();
    virtual ~PersEx01EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
