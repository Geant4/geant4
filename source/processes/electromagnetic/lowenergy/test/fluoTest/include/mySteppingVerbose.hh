//---------------------------------------------------------------
//
// mySteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class mySteppingVerbose;

#ifndef mySteppingVerbose_h
#define mySteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class mySteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  mySteppingVerbose();
 ~mySteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif

