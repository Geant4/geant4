#ifndef MLRunManager_h
#define MLRunManager_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4RunManager.hh"
#include "globals.hh"

#include "G4Timer.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLRunManager: public G4RunManager {

public:
  MLRunManager ();

  virtual ~MLRunManager ();

public:

  G4double GetTimeUsed ();


public:
  // to force the re-initialization of physicsList() 
  //  inline void CutOffHasNotBeenChanged (G4bool bval)
  //{ cutoffInitialized = bval; };

  inline void PhysicsListHasNotBeenChanged (G4bool bval) 
  { physicsInitialized = bval; };

private:

  G4Timer* timer;
  G4double timeused;

};
////////////////////////////////////////////////////////////////////////////////
//
#endif
