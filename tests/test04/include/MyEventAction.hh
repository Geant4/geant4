
#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction();
    virtual ~MyEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event* );

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
