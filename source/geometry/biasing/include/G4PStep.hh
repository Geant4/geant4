#ifndef G4PStep_hh
#define G4PStep_hh G4PStep_hh 

#include "G4PTouchableKey.hh"

class G4PStep {
public:
  G4PStep(const G4PTouchableKey &preKey, const G4PTouchableKey &postKey):
    fPreTouchableKey(preKey), 
    fPostTouchableKey(postKey){}
  ~G4PStep(){}
  G4PTouchableKey fPreTouchableKey;
  G4PTouchableKey fPostTouchableKey;  
  G4bool fCrossBoundary;
};


#endif
