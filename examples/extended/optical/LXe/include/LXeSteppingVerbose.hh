
class LXeSteppingVerbose;

#ifndef LXeSteppingVerbose_h
#define LXeSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"


class LXeSteppingVerbose : public G4SteppingVerbose
{
 public:   

   LXeSteppingVerbose();
  ~LXeSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
