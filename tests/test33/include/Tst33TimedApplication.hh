#ifndef Tst33TimedApplication_hh
#define Tst33TimedApplication_hh Tst33TimedApplication_hh

#include "Tst33VApplication.hh"
#include "globals.hh"

class Tst33TimedEventAction;


class Tst33TimedApplication : public Tst33VApplication {
public:
  explicit Tst33TimedApplication(G4int time);
  virtual ~Tst33TimedApplication();
  
  virtual G4UserRunAction *CreateRunAction();
  virtual Tst33VEventAction *CreateEventAction();
private:
  G4int fTime;
};

#endif
