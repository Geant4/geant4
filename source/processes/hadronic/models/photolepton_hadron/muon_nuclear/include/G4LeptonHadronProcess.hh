//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
      virtual ~G4LeptonHadronProcess();

      virtual G4VParticleChange *PostStepDoIt( const G4Track &leptonTrack,
					       const G4Step &aStep );

    protected:

      G4LeptonHadronInteractionModel* theInteractionModel;
      G4Nucleus targetNucleus;

  };

#endif
