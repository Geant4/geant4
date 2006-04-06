// -------------------------------------------------------------------
// $Id: MicrobeamSteppingVerbose.hh,v 1.1 2006-04-06 15:32:44 sincerti Exp $
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
