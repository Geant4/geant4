//---------------------------------------------------------------
//
// fluoTestSteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class fluoTestSteppingVerbose;

#ifndef fluoTestSteppingVerbose_h
#define fluoTestSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class fluoTestSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  fluoTestSteppingVerbose();
 ~fluoTestSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif

