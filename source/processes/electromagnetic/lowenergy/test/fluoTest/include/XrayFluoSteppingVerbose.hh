//---------------------------------------------------------------
//
// XrayFluoSteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class XrayFluoSteppingVerbose;

#ifndef XrayFluoSteppingVerbose_h
#define XrayFluoSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class XrayFluoSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  XrayFluoSteppingVerbose();
 ~XrayFluoSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif

