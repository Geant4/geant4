#ifndef MyVisManager_h
#define MyVisManager_h 1

#include "G4VisManager.hh"

class MyVisManager: public G4VisManager {

public:

  MyVisManager();

private:

  void RegisterGraphicsSystems ();

};

#endif
