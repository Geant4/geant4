
#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction();
    ~MyEventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
