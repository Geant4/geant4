#include "G4ParallelImportanceProcess.hh"
#include "G4VImportanceSampler.hh"
#include "g4std/strstream"

G4ParallelImportanceProcess::
G4ParallelImportanceProcess(const G4VImportanceSampler &aImportanceSampler,
		    G4VPGeoDriver &pgeodriver,
		    G4VParallelStepper &aStepper, 
		    const G4String &aName):
  G4ParallelTransport(pgeodriver, aStepper, aName),
  fParticleChange(G4ParallelTransport::fParticleChange),
  fImportanceSampler(aImportanceSampler)
{}


G4VParticleChange *G4ParallelImportanceProcess::
PostStepDoIt(const G4Track& aTrack, const G4Step &aStep){
  if (aTrack.GetTrackStatus()==fStopAndKill) {
    G4cout << "G4ParallelImportanceProcess::PostStepDoIt StopAndKill" << G4endl;
  }
  G4ParallelTransport::PostStepDoIt(aTrack, aStep);

  // get new weight and number of clones
  G4Nsplit_Weight nw(fImportanceSampler.Sample(aTrack.GetWeight()));

  fImportancePostStepDoIt.DoIt(aTrack, fParticleChange, nw);
  return fParticleChange;
}
  
