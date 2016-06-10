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
// $Id: G4ITSafetyHelper.cc 72309 2013-07-15 15:52:17Z gcosmo $
// GEANT4 tag $ Name:  $
// 
// class G4ITSafetyHelper Implementation
//
// Original author: John Apostolakis, 2006
//
// --------------------------------------------------------------------

#include "G4ITSafetyHelper.hh"
#include "G4ITPathFinder.hh"
#include "G4ITTransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4PathFinder.hh"
#include "globals.hh"

G4ITSafetyHelper::G4ITSafetyHelper() :
    G4TrackStateDependent<G4ITSafetyHelper>(), fUseParallelGeometries(false), // By default, one geometry only
    fFirstCall(true), fVerbose(0)
// fRecomputeFactor(0.0)
{
  fpPathFinder = 0; //  Cannot initialise this yet - a loop results

  // Initialization of the Navigator pointer is postponed, and must
  // be undertaken by another class calling InitialiseHelper()
  //
  fpMassNavigator = 0;
  fMassNavigatorId = -1;
}

void G4ITSafetyHelper::InitialiseNavigator()
{
  fpPathFinder = G4PathFinder::GetInstance();

  G4ITTransportationManager* pTransportMgr =
      G4ITTransportationManager::GetTransportationManager();

  fpMassNavigator = pTransportMgr->GetNavigatorForTracking();

  if(fpMassNavigator == 0) abort();

  // Check
  //
  G4VPhysicalVolume* worldPV = fpMassNavigator->GetWorldVolume();
  if (worldPV == 0)
  {
    G4Exception("G4ITSafetyHelper::InitialiseNavigator",
                "InvalidNavigatorWorld", FatalException,
                "Found that existing tracking Navigator has NULL world");
  }

  // fMassNavigatorId = pTransportMgr->ActivateNavigator( fpMassNavigator );
}

void G4ITSafetyHelper::InitialiseHelper()
{
  NewTrackState();
  if (fFirstCall)
  {
    InitialiseNavigator();
  }
  fFirstCall = false;
}

G4ITSafetyHelper::~G4ITSafetyHelper()
{
}

G4double G4ITSafetyHelper::CheckNextStep(const G4ThreeVector &position,
                                         const G4ThreeVector &direction,
                                         const G4double currentMaxStep,
                                         G4double& newSafety)
{
  // Distance in the Mass geometry
  //
  G4double linstep = fpMassNavigator->CheckNextStep(position, direction,
                                                    currentMaxStep, newSafety);
  fpTrackState->fLastSafetyPosition = position;
  fpTrackState->fLastSafety = newSafety;

  // TODO: Can replace this with a call to PathFinder
  //        giving id of Mass Geometry --> this avoid doing the work twice

  return linstep;
}

G4double G4ITSafetyHelper::ComputeSafety(const G4ThreeVector& position,
                                         G4double maxLength)
{
  G4double newSafety;

  // Only recompute (calling Navigator/PathFinder) if 'position'
  // is  *not* the safety location and has moved 'significantly'
  //
  G4double moveLengthSq = (position - fpTrackState->fLastSafetyPosition).mag2();
  if ((moveLengthSq > 0.0))
  {
    if (!fUseParallelGeometries)
    {
      // Safety for mass geometry
      newSafety = fpMassNavigator->ComputeSafety(position, maxLength, true);
    }
    else
    {
      // Safety for all geometries
      newSafety = fpPathFinder->ComputeSafety(position);
    }

    // We can only store a 'true' safety - one that was not restricted by maxLength
    if (newSafety < maxLength)
    {
      fpTrackState->fLastSafety = newSafety;
      fpTrackState->fLastSafetyPosition = position;
    }
  }
  else
  {
    // return last value if position is not (significantly) changed
    //
    // G4double moveLength = 0;
    // if( moveLengthSq > 0.0 ) { moveLength= std::sqrt(moveLengthSq); }
    newSafety = fpTrackState->fLastSafety; // -moveLength;
  }
  return newSafety;
}

void G4ITSafetyHelper::ReLocateWithinVolume(const G4ThreeVector &newPosition)
{
#ifdef G4VERBOSE
  if (fVerbose > 0)
  {
    // There is an opportunity - and need - to check whether the proposed move is safe
    G4ThreeVector moveVec = newPosition - fpTrackState->fLastSafetyPosition;
    if (moveVec.mag2() > sqr(fpTrackState->fLastSafety))
    {
      // A problem exists - we are proposing to move outside 'Safety Sphere'
      G4ExceptionDescription ed;
      ed << " Safety Sphere:  Radius = " << fpTrackState->fLastSafety;
      ed << " Center   = " << fpTrackState->fLastSafetyPosition << G4endl;
      ed << " New Location :  Move   = " << moveVec.mag2();
      ed << " Position = " << newPosition << G4endl;
      G4Exception("G4ITSafetyHelper::ReLocateWithinVolume", "GeomNav999",
                  JustWarning,
                  "Unsafe Move> Asked to relocate beyond 'Safety sphere'.");
    }
  }
#endif

  if (!fUseParallelGeometries)
  {
    fpMassNavigator->LocateGlobalPointWithinVolume(newPosition);
  }
  else
  {
    fpPathFinder->ReLocate(newPosition);
  }
}

void G4ITSafetyHelper::Locate(const G4ThreeVector& newPosition,
                              const G4ThreeVector& newDirection)
{
  if (!fUseParallelGeometries)
  {
    fpMassNavigator->LocateGlobalPointAndSetup(newPosition, &newDirection, true,
                                               false);
  }
  else
  {
    fpPathFinder->Locate(newPosition, newDirection);
  }
}
