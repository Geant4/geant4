#ifndef B01VisApplication_hh
#define B01VisApplication_hh B01VisApplication_hh

#include "B01VApplication.hh"
#include "B01VisRun.hh"
class G4UIsession;

class B01VisApplication : public B01VApplication {
public:
  B01VisApplication();
  virtual ~B01VisApplication();
  virtual void RunSimulation(B01VSimulation *sim);
  G4UIsession *CreateSession();
private:
  B01VisRun fVisRun;
};

#endif
