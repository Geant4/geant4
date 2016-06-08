
#ifndef ExN04EventAction_h
#define ExN04EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ExN04EventAction : public G4UserEventAction
{
  public:
    ExN04EventAction();
    ~ExN04EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    G4int trackerCollID;
    G4int calorimeterCollID;
    G4int muonCollID;
};

#endif

    
