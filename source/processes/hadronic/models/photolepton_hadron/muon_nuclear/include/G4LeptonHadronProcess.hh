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
// G4LeptonHadronProcess.hh
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#ifndef G4LeptonHadronProcess_h
#define G4LeptonHadronProcess_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChange.hh"
#include "G4LeptonHadronInteractionModel.hh"


  class G4LeptonHadronProcess : public G4VDiscreteProcess
  {
    public:

      G4LeptonHadronProcess(const G4String &processName = "LeptonHadron");
      ~G4LeptonHadronProcess();

      G4VParticleChange *PostStepDoIt( const G4Track &leptonTrack,
                                       const G4Step &aStep )
      {
        return G4LeptonHadronProcess::GeneralPostStepDoIt(leptonTrack, aStep);
      }

      G4VParticleChange *GeneralPostStepDoIt( const G4Track &leptonTrack,
                                              const G4Step &aStep );

      virtual G4LeptonHadronInteractionModel *chooseInteractionModel() = 0;


    protected:

      G4LeptonHadronInteractionModel* theInteractionModel;
      G4Nucleus targetNucleus;

  };

#endif
