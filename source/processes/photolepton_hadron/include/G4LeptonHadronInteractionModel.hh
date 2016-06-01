// G4LeptonHadronInteractionModel.hh
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#ifndef G4LeptonHadronInteractionModel_h
#define G4LeptonHadronInteractionModel_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4ParticleChange.hh"


  class G4LeptonHadronInteractionModel
  {
    public:

      G4LeptonHadronInteractionModel();
      ~G4LeptonHadronInteractionModel();

      virtual void makePhysicsVector() = 0;
      virtual G4VParticleChange* applyInteractionModel
        (const G4Track &leptonTrack, G4Nucleus &targetNucleus) = 0;
      virtual G4double computeMicroscopicCrossSection
                             (const G4Track &leptonTrack) = 0;

      inline G4int operator==(
        const G4LeptonHadronInteractionModel &right) const
      { return this == &right; }

      inline G4int operator!=(
        const G4LeptonHadronInteractionModel &right) const
      { return this != &right; }


    protected:

      G4ParticleChange aParticleChange;


    private:

      inline G4LeptonHadronInteractionModel
        (const G4LeptonHadronInteractionModel &right) { }

      inline G4LeptonHadronInteractionModel & operator=(
        const G4LeptonHadronInteractionModel &right)
      { return *this; }

  };

#endif
