#ifndef Tst22EventAction_h
#define Tst22EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class Tst22EventAction : public G4UserEventAction
{
  public:
    Tst22EventAction(){evnum=0;}
    ~Tst22EventAction(){}

  public:
    virtual void BeginOfEventAction(const G4Event*)
    {
      G4cout <<"Event number "<<evnum++<<G4endl;
    }
    virtual void EndOfEventAction(const G4Event*){ }
 
  private:
    int evnum;
};

#endif

    
