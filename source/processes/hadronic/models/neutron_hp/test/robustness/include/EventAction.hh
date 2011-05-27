#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <iostream>
#include "G4Event.hh"

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction(); 

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    
  private:
    int nEvent;
};

#endif

    
