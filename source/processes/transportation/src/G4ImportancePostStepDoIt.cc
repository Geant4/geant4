#include "G4ImportancePostStepDoIt.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"

#include "G4VImportanceSampler.hh"
#include "G4Nsplit_Weight.hh"

#include <strstream>

void G4ImportancePostStepDoIt::DoIt(const G4Track& aTrack, 
				    G4ParticleChange *aParticleChange,
				    const G4Nsplit_Weight nw){
  
  // evaluate results from sampler
  if (nw.fN>1) {
    // split track 
    Split(aTrack, nw, aParticleChange);
  }
  else if (nw.fN==1) {
    // don't split, but weight may be changed ! 
    aParticleChange->SetWeightChange(nw.fW);
  }
  else if (nw.fN==0) {
    // kill track
    aParticleChange->SetStatusChange(fStopAndKill);
  }
  else {
    // wrong answer
    G4std::ostrstream os;
    os << "G4ImportancePostStepDoIt::DoIt: sampler returned nw = " 
       << nw << '\0' << G4endl;
    G4Exception(os.str());
  }
}






void G4ImportancePostStepDoIt::Split(const G4Track &aTrack,
				     const G4Nsplit_Weight &nw,
				     G4ParticleChange *aParticleChange) {
  aParticleChange->SetWeightChange(nw.fW);
  aParticleChange->SetNumberOfSecondaries(nw.fN-1);
  
  for (G4int i=1;i<nw.fN;i++) {
    G4Track *ptrack = new G4Track(aTrack);
    
    //    ptrack->SetCreatorProcess(aTrack.GetCreatorProcess());
    ptrack->SetWeight(nw.fW);
    
    if (ptrack->GetMomentumDirection() != aTrack.GetMomentumDirection()) {
      G4Exception("ERROR:G4ImportancePostStepDoIt::Split: (ptrack->GetMomentumDirection() != aTrack.GetMomentumDirection()");
    }
    
    aParticleChange->AddSecondary(ptrack);
  }
  return;
}  
