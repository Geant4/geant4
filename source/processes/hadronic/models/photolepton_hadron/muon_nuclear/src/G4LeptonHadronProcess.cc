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
// G4LeptonHadronProcess.cc
//
//     M.Takahata (Makoto.Takahata@cern.ch)

#include "G4LeptonHadronProcess.hh"

//-----------------------------------------------------------------------------
  G4LeptonHadronProcess::G4LeptonHadronProcess( const G4String &processName )
//-----------------------------------------------------------------------------
    : G4VDiscreteProcess( processName )
  {
  }


//-----------------------------------------------------------------------------
  G4LeptonHadronProcess::~G4LeptonHadronProcess()
//-----------------------------------------------------------------------------
  {
  }


//-----------------------------------------------------------------------------
  G4VParticleChange* 
  G4LeptonHadronProcess::GeneralPostStepDoIt( const G4Track &leptonTrack,
                                              const G4Step & )
//-----------------------------------------------------------------------------
  {
    targetNucleus.ChooseParameters(leptonTrack.GetMaterial());
    G4VParticleChange *result 
      = theInteractionModel->applyInteractionModel(leptonTrack, targetNucleus);

    ResetNumberOfInteractionLengthLeft();

    return result;
  }
