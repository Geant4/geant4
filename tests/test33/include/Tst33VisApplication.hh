#ifndef Tst33VisApplication_hh
#define Tst33VisApplication_hh Tst33VisApplication_hh

#include "Tst33VApplication.hh"
#include "Tst33VisManager.hh"




class Tst33VisApplication : public Tst33VApplication {
public:
  Tst33VisApplication();
  virtual ~Tst33VisApplication();

  virtual G4UserRunAction *CreateRunAction();
  virtual Tst33VEventAction *CreateEventAction();
  

private:
  Tst33VisManager fVisManager;
};

#endif
