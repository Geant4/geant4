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
