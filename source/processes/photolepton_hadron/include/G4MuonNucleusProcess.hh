// G4MuonNucleusProcess.hh
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#ifndef G4MuonNucleusProcess_h
#define G4MuonNucleusProcess_h 1

#include "globals.hh"
#include "G4LeptonHadronProcess.hh"
#include "G4MuonNucleusInteractionModel.hh"


  class G4MuonNucleusProcess : public G4LeptonHadronProcess
  {
    public:

      G4MuonNucleusProcess(const G4String& processName ="MuonNucleus");
      ~G4MuonNucleusProcess();


      G4double GetMeanFreePath(const G4Track &muonTrack,
                               G4double previousStepSize,
                               G4ForceCondition *condition);

      G4LeptonHadronInteractionModel *chooseInteractionModel();

      G4VParticleChange *PostStepDoIt(const G4Track &muonTrack,
                                      const G4Step &aStep)
      { 
        return G4LeptonHadronProcess::GeneralPostStepDoIt(muonTrack, aStep);
      }


    private:
      G4MuonNucleusProcess(const G4MuonNucleusProcess &right);

  };

#endif
