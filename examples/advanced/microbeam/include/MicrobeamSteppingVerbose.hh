// -------------------------------------------------------------------
// $Id: MicrobeamSteppingVerbose.hh,v 1.3 2006-06-01 22:25:19 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamSteppingVerbose_h
#define MicrobeamSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class MicrobeamSteppingVerbose : public G4SteppingVerbose {

public:   
  
  MicrobeamSteppingVerbose();
  ~MicrobeamSteppingVerbose();
  
  void StepInfo();
  void TrackingStarted();

};

#endif
