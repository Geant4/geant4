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
//
// -- Bogus -- BaBar Object-Oriented Geant-based Unified Simulation
//
// BgsLooperDeath
//
// Description:
//   This is a simple GEANT4 process that destroys any particle below
//   the specified total kinetic energy that reverse direction (loops 180
//   degress) in the x/y plane.
//
//   Based on BgsChargedLowEnergyDeath
//
// Author List:
//   David Williams
//
// Modification History:
//
//-----------------------------------------------------------------------------
#ifndef NTSTlooperDeath_hh
#define NTSTLooperDeath_hh 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VProcess.hh"

class NTSTLooperDeath : public G4VProcess
{
public:
  
  NTSTLooperDeath( G4double theMinMomentum=5*CLHEP::MeV,
		   const char* name="NTSTLoopDeath",
		   G4ProcessType type=fUserDefined );

  ~NTSTLooperDeath();
	
	
  //
  // Derived methods
  //
  virtual G4double
  PostStepGetPhysicalInteractionLength( const G4Track& track,
					G4double   previousStepSize,
					G4ForceCondition* condition );

  virtual G4VParticleChange* PostStepDoIt( const G4Track &track,
					   const G4Step &step );

  virtual G4double
  AlongStepGetPhysicalInteractionLength( const G4Track&,
					 G4double  , // previousStepSize,
					 G4double  , // currentMinimumStep,
					 G4double& , // currentSafety,
					 G4GPILSelection* ) // selection )
    { return -1.0; }

  virtual G4VParticleChange* AlongStepDoIt( const G4Track & , // track,
					    const G4Step & ) // step
    { return 0; }
 
  virtual G4double
  AtRestGetPhysicalInteractionLength( const G4Track &,          // track,
				      G4ForceCondition * )      // force )
    { return -1.0; }

  virtual G4VParticleChange*
  AtRestDoIt( const G4Track &, const G4Step & )     // track, step )
    { return 0; }


  virtual G4bool IsApplicable( const G4ParticleDefinition &particle )
    { return (particle.GetPDGCharge() != 0); }
	
  //
  // Accessors
  //
  void SetMinMomentum( G4double theMinMomentum )
    { minMomentum = theMinMomentum; }
  G4double GetMinMomentum() const
    { return minMomentum; }

protected:
  
  G4double minMomentum; 
};

#endif
