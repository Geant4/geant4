
#ifndef ExampleVisManager_h
#define ExampleVisManager_h 1

#include "G4VisManager.hh"


class ExampleVisManager: public G4VisManager {

public:

  ExampleVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif
