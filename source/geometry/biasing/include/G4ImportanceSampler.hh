#ifndef G4ImportanceSampler_hh
#define G4ImportanceSampler_hh G4ImportanceSampler_hh

#include "G4VImportanceSampler.hh"

class G4VImportanceAlgorithm;
class G4VParallelStepper;
class G4ImportanceFinder;
class G4VIStore;

class G4ImportanceSampler: public G4VImportanceSampler {
public:
  G4ImportanceSampler(const G4VImportanceAlgorithm &aIalg,
		      const G4VParallelStepper &astepper,
		      const G4VIStore &istore);
  ~G4ImportanceSampler();
  G4Nsplit_Weight Sample(G4double w) const; 
private:
  G4ImportanceSampler(const G4ImportanceSampler &);
  G4ImportanceSampler &operator=(const G4ImportanceSampler &);

  const G4VImportanceAlgorithm &fIalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4ImportanceFinder &fIfinder;
};


#endif
