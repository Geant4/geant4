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
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History:
//      17 August 2004  P. Gumplinger, T. MacPhail
// ------------------------------------------------------------
//
#ifndef G4DecayWithSpin_h
#define G4DecayWithSpin_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"

#include "G4Decay.hh"

class G4VExtDecayer;

class G4DecayWithSpin : public G4Decay
{
  public:
    //  Constructors
    G4DecayWithSpin(const G4String& processName ="DecayWithSpin");

    //  Destructor
    virtual ~G4DecayWithSpin();

  protected: // With Description
    virtual G4VParticleChange* DecayIt(
                             const G4Track& aTrack,
                             const G4Step&  aStep
                            );
    // The DecayIt() method returns by pointer a particle-change object,
    // which has information of daughter particles.

  private:
  G4ThreeVector Spin_Precession(const G4Step& aStep,
                                G4ThreeVector B, G4double deltatime );

};

#endif
