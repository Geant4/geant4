#ifndef B01TimedApplication_hh
#define B01TimedApplication_hh B01TimedApplication_hh

#include "B01VApplication.hh"

#include "B01TimedRun.hh"

class B01TimedApplication : public B01VApplication {
public:
  explicit B01TimedApplication(G4int time);
  virtual ~B01TimedApplication();
  virtual void RunSimulation(B01VSimulation *sim);
private:
  B01TimedRun fTimedRun;
  
};

#endif
