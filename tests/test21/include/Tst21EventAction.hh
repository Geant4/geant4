#ifndef Tst21EventAction_h
#define Tst21EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class Tst21EventAction : public G4UserEventAction
{
  public:
    Tst21EventAction(){evnum=0;}
    ~Tst21EventAction(){}

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

    
