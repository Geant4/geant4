#ifndef HECEventAction_h
#define HECEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class HECEventAction : public G4UserEventAction
{
  public:
    HECEventAction();
    ~HECEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    void BadEvent(){GoodEv = false;}
 
  private:
    static int evnum;
    static G4bool GoodEv;
    static G4bool store;
    G4int moveCollID;
    G4int frontCollID;
    G4int larCollID;
    G4int cuCollID;
    G4int leakCollID;
};

#endif

    
