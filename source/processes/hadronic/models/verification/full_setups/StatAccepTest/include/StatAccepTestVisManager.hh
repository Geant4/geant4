#ifndef StatAccepTestVisManager_h
#define StatAccepTestVisManager_h 1

#include "G4VisManager.hh"

class StatAccepTestVisManager: public G4VisManager {

public:

  StatAccepTestVisManager();

private:

  void RegisterGraphicsSystems ();

};

#endif
