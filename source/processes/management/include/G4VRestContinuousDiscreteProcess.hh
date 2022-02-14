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
// G4VRestContinuousDiscreteProcess
//
// Class description:
//
// Abstract class which defines the public behavior of
// continuous and discrete physics interactions.

// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 8 January 1997, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------
#ifndef G4VRestContinuousDiscreteProcess_hh
#define G4VRestContinuousDiscreteProcess_hh 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VProcess.hh"

class G4VRestContinuousDiscreteProcess : public G4VProcess 
{
  public:     

    G4VRestContinuousDiscreteProcess(const G4String& ,
                                     G4ProcessType aType = fNotDefined);
    G4VRestContinuousDiscreteProcess(G4VRestContinuousDiscreteProcess &);

    virtual ~G4VRestContinuousDiscreteProcess();

    G4VRestContinuousDiscreteProcess& operator=(const G4VRestContinuousDiscreteProcess&) = delete;

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

    virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );

    virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4double currentMinimumStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection
                            );

    virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& ,
                             const G4Step& 
                            );
 
    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4ForceCondition* 
                            );

    virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            );

  protected:

    virtual G4double GetMeanLifeTime( const G4Track& aTrack,
                                      G4ForceCondition* condition ) = 0;
      // Calculates the mean life-time (i.e. for decays) of the
      // particle at rest due to the occurrence of the given process,
      // or converts the probability of interaction (i.e. for
      // annihilation) into the life-time of the particle for the
      // occurrence of the given process

    virtual G4double GetContinuousStepLimit( const G4Track& aTrack,
                                             G4double  previousStepSize,
                                             G4double  currentMinimumStep,
                                             G4double& currentSafety ) = 0;

    virtual G4double GetMeanFreePath( const G4Track& aTrack,
                                      G4double   previousStepSize,
                                      G4ForceCondition* condition ) = 0;
      // Calculates from the macroscopic cross-section a mean
      // free path, the value is returned in units of distance

    inline void SetGPILSelection(G4GPILSelection selection)
      { valueGPILSelection = selection; }

    inline G4GPILSelection GetGPILSelection() const
      { return valueGPILSelection; }

  private:

    G4VRestContinuousDiscreteProcess();
      // Hidden default constructor

    G4GPILSelection valueGPILSelection = CandidateForSelection;
      // This is the returned value of G4GPILSelection in 
      // the arguments of AlongStepGPIL()
};

#endif
