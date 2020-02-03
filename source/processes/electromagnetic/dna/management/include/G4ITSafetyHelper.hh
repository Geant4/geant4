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
//
//
// class G4ITSafetyHelper
//
// Class description:
//
// G4ITSafetyHelper is a duplicated version of G4SafetyHelper
// This class is a helper for physics processes which require 
// knowledge of the safety, and the step size for the 'mass' geometry

// First version:  J. Apostolakis,  July 5th, 2006
// Modified:
//  10.04.07 V.Ivanchenko  Use unique G4ITSafetyHelper
// --------------------------------------------------------------------

#ifndef G4ITSAFETYHELPER_HH
#define G4ITSAFETYHELPER_HH 1

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4ITNavigator.hh"
#include "G4TrackState.hh"

class G4PathFinder;

class G4ITSafetyHelper : public G4TrackStateDependent<G4ITSafetyHelper>
{
public:
  // with description
  G4ITSafetyHelper();
  ~G4ITSafetyHelper();
  //
  // Constructor and destructor

  G4double CheckNextStep(const G4ThreeVector& position,
                         const G4ThreeVector& direction,
                         const G4double currentMaxStep,
                         G4double& newSafety);
  //
  // Return linear step for mass geometry

  G4double ComputeSafety(const G4ThreeVector& pGlobalPoint,
                         G4double maxRadius = DBL_MAX); // Radius of interest
  //
  // Return safety for all geometries.
  //
  //  The 2nd argument is the radius of your interest (e.g. maximum displacement )
  //    Giving this you can reduce the average computational cost.
  //  If the second argument is not given, this is the real isotropic safety

  void Locate(const G4ThreeVector& pGlobalPoint,
              const G4ThreeVector& direction);
  //
  // Locate the point for all geometries

  void ReLocateWithinVolume(const G4ThreeVector& pGlobalPoint );
  //
  // Relocate the point in the volume of interest

  inline void EnableParallelNavigation(G4bool parallel);
  //
  //  To have parallel worlds considered, must be true.
  //  Alternative is to use single (mass) Navigator directly

  void InitialiseNavigator();
  //
  // Check for new navigator for tracking, and reinitialise pointer

  G4int SetVerboseLevel( G4int lev )
  { G4int oldlv= fVerbose; fVerbose= lev; return oldlv;}

  inline G4VPhysicalVolume* GetWorldVolume();
  inline void SetCurrentSafety(G4double val, const G4ThreeVector& pos);

public: // without description

  void InitialiseHelper();

private:

  G4PathFinder* fpPathFinder;
  G4ITNavigator* fpMassNavigator;
  G4int fMassNavigatorId;

  G4bool fUseParallelGeometries;
  // Flag whether to use PathFinder or single (mass) Navigator directly
  G4bool fFirstCall;
  // Flag of first call
  G4int fVerbose;
  // Whether to print warning in case of move outside safety

public:
  // State used during tracking -- for optimisation
  class State
  {
    friend class G4ITSafetyHelper;
    G4ThreeVector fLastSafetyPosition;
    G4double fLastSafety;

  public:
    State() :
    fLastSafetyPosition(0.0,0.0,0.0),
    fLastSafety(0.0)
    {}

    virtual ~State()
    {}
  };

  // const G4double  fRecomputeFactor;
  // parameter for further optimisation:
  // if ( move < fact*safety )  do fast recomputation of safety
  // End State (tracking)
};

RegisterTrackState(G4ITSafetyHelper, State)

// Inline definitions

inline
void G4ITSafetyHelper::EnableParallelNavigation(G4bool parallel)
{
  fUseParallelGeometries = parallel;
}

inline G4VPhysicalVolume* G4ITSafetyHelper::GetWorldVolume()
{
  return fpMassNavigator->GetWorldVolume();
}

inline
void G4ITSafetyHelper::SetCurrentSafety(G4double val, const G4ThreeVector& pos)
{
  fpTrackState->fLastSafety = val;
  fpTrackState->fLastSafetyPosition = pos;
}

#endif
