#ifndef G4ParallelStepper_hh
#define G4ParallelStepper_hh G4ParallelStepper_hh

#include "G4VParallelStepper.hh"
#include "G4PStep.hh"

class G4ParallelStepper : public G4VParallelStepper {
public:
  G4ParallelStepper();
  ~G4ParallelStepper();
  G4ParallelStepper(const G4ParallelStepper &);
  G4ParallelStepper &operator=(const G4ParallelStepper &);

  G4PStep GetPStep() const {return *fPStep;}
  void Init(const G4PTouchableKey &aptk);
  void Update(const G4PTouchableKey &aptk);
  void UnSetCrossBoundary();
private:
  G4PStep *fPStep;
  void Error(const G4String &m) {
    G4cout << "ERROE: in G4ParallelStepper::" <<  m << G4endl;
    exit(1);
  }
};


#endif



