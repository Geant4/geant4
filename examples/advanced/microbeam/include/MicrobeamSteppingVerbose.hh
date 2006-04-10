// -------------------------------------------------------------------
// $Id: MicrobeamSteppingVerbose.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
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
