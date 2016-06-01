
#ifndef ExE02EventAction_h
#define ExE02EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ExE02EventAction : public G4UserEventAction
{
  public:
    ExE02EventAction();
    ~ExE02EventAction();

  public:
    void BeginOfEventAction();
    void EndOfEventAction();

  private:
    G4int colID1;
    G4int colID2;
};

#endif

    
