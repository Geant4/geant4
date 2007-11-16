#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class EventAction : public G4UserEventAction {

 public:
   EventAction();
   ~EventAction();

   void BeginOfEventAction(const G4Event*);
   void EndOfEventAction(const G4Event*);
};

#endif // EVENTACTION_HH
