#include "Tst33TimedApplication.hh"
#include "G4RunManager.hh"
#include "Tst33TimedEventAction.hh"



Tst33TimedApplication::Tst33TimedApplication(G4int time)
  :
  fTime(time)
{}

Tst33TimedApplication::~Tst33TimedApplication()
{}

Tst33VEventAction *Tst33TimedApplication::CreateEventAction() {
  return new Tst33TimedEventAction(fTime);
}



G4UserRunAction *Tst33TimedApplication::CreateRunAction(){
  return 0;
}

