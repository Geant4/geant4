#include "EventAction.hh"
    EventAction::EventAction() { nEvent = 0;}
    EventAction::~EventAction(){}

    void   EventAction::BeginOfEventAction(const G4Event*)
    {
      nEvent++;
      cout << "Processing event "<<nEvent<<endl;
    }
