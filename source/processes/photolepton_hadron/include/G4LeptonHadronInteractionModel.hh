//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
      virtual ~G4LeptonHadronInteractionModel();

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
