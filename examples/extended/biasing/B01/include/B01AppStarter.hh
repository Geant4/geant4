#ifndef B01AppStarter_hh
#define B01AppStarter_hh B01AppStarter_hh

#include "globals.hh"
#include "B01AppStarterMessenger.hh"

class B01VApplication;
class B01VSimulation;


class B01AppStarter {
public:
  B01AppStarter();
  ~B01AppStarter();
  void SetSimulation(B01VSimulation *sim);
  void CreateVisApplication();
  void CreateTimedApplication(G4int timed);
  void Run();
private:
  B01AppStarter(const B01AppStarter &);
  B01AppStarter &operator=(const B01AppStarter &);

  B01AppStarterMessenger fMessenger;
  B01VSimulation *fSim;
  B01VApplication *fApp;

};

#endif
