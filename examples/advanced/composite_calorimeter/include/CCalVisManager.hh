///////////////////////////////////////////////////////////////////////////////
// File: CCalVisManager.hh
// Description: Example Visualization Manager implementing virtual function
//              RegisterGraphicsSystems.  Exploits C-pre-processor variables
//              G4VIS_USE_DAWN, etc., which are set by the GNUmakefiles if
//              environment variables of the same name are set.
//
//              So all you have to do is set environment variables and compile
//              and instantiate this in your main().
//
//              Alternatively, you can implement an empty function here and 
//              just register the systems you want in your main(), e.g.:
//              G4VisManager* myVisManager = new CCalVisManager;
//              myVisManager -> RegisterGraphicsSystem (new MyGraphicsSystem);
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalEXAMPLEVISMANAGER_HH
#define CCalEXAMPLEVISMANAGER_HH

#ifdef G4VIS_USE
#include "G4VisManager.hh"

class CCalVisManager: public G4VisManager {

public:

  CCalVisManager (G4int verboseLevel = 0);
  // Controls initial verbose level of VisManager and VisMessenger.
  // Can be changed by /vis/set/verbose.

private:

  virtual void RegisterGraphicsSystems ();

};

#endif
#endif
