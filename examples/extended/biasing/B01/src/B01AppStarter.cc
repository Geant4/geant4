#include "B01AppStarter.hh"
#include "B01SimulationFactory.hh"
#include "B01VSimulation.hh"
#include "B01VisApplication.hh"
#include "B01TimedApplication.hh"

B01AppStarter::B01AppStarter()
  : 
  fMessenger(this),
  fSim(0),
  fApp(0)
{}
B01AppStarter::~B01AppStarter()
{
  if (fApp) {
    delete fApp;
  }
}

void B01AppStarter::SetSimulation(B01VSimulation *sim){
  fSim = sim;
}

void B01AppStarter::CreateVisApplication(){
  fApp = new B01VisApplication;
  if (!fApp) {
    G4std::G4Exception("B01AppStarter::CreateVisApplication: new failed to create B01VisApplication!");
  }
}

void B01AppStarter::CreateTimedApplication(G4int time) {
  fApp = new B01TimedApplication(time);
  if (!fApp) {
    G4std::G4Exception("B01AppStarter::CreateTimedApplication: new failed to create B01TimedApplication!");
  }
}

void B01AppStarter::Run(){
  if (!fSim || !fApp) {
    if (!fSim) {
      G4std::G4cout << "B01AppStarter::Run(): Error: create a simulation first!" << G4endl;    
    }
    if (!fApp) {
      G4std::G4cout << "B01AppStarter::Run(): Error: create an application first!" << G4endl;
    }    
  }
  else {
    fApp->RunSimulation(fSim);
  }
}
