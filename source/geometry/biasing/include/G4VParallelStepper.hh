#ifndef G4VParallelStepper_hh
#define G4VParallelStepper_hh G4VParallelStepper_hh

#include "G4PStep.hh"

class G4VParallelStepper {
public:
  virtual ~G4VParallelStepper(){}
  virtual G4PStep GetPStep() const = 0;
  virtual void Init(const G4PTouchableKey &aptk) = 0;
  virtual void Update(const G4PTouchableKey &aptk) = 0;
  virtual void UnSetCrossBoundary() = 0;
};


#endif
