#include "G4ImportanceSampler.hh"
#include "G4ImportanceFinder.hh"
#include "G4VParallelStepper.hh"
#include "G4VImportanceAlgorithm.hh"

G4ImportanceSampler::
G4ImportanceSampler(const G4VImportanceAlgorithm &aIalg,
		    const G4VParallelStepper &astepper,
		    const G4VIStore &istore):
  fIalgorithm(aIalg),
  fPStepper(astepper),
  fIfinder(*(new G4ImportanceFinder(istore))){}

G4ImportanceSampler::~G4ImportanceSampler(){
  delete &fIfinder;
}
  
G4Nsplit_Weight
G4ImportanceSampler::
Sample(G4double w) const {
  G4PStep pstep = fPStepper.GetPStep();
  return fIalgorithm.
    Calculate(fIfinder.GetIPre_over_IPost(pstep.fPreTouchableKey,
					  pstep.fPostTouchableKey), 
	      w);
}

