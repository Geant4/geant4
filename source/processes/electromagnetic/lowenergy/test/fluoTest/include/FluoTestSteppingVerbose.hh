//---------------------------------------------------------------
//
// FluoTestSteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class FluoTestSteppingVerbose;

#ifndef FluoTestSteppingVerbose_h
#define FluoTestSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class FluoTestSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  FluoTestSteppingVerbose();
 ~FluoTestSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif

