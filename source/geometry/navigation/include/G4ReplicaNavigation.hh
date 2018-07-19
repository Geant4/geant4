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
// $Id: G4ReplicaNavigation.hh 96458 2016-04-15 10:15:24Z gcosmo $
//
// 
// class G4ReplicaNavigation
//
// Class description:
//
// Utility for navigation in volumes containing a single G4PVParameterised
// volume for which voxels for the replicated volumes have been constructed.
// [Voxels MUST be along one axis only: NOT refined]

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef G4REPLICANAVIGATION_HH
#define G4REPLICANAVIGATION_HH

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4Types.hh"
#include "G4NavigationHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4BlockingList.hh"

// Required for voxel handling
//
#include "G4SmartVoxelHeader.hh"

class G4VSolid;

struct G4ExitNormal
{
   // Bucket to hold value of Normal (3-vector), 
   // bools for calculated and leave-behind or 'validConvex',
   // and exiting side.
   
   enum  ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPX,kMX,kPY,kMY,kPZ,kMZ,kMother};
     // Identity of 'Side' of Replicas. Used by DistanceToOut methods.

   G4ThreeVector exitNormal;
   G4bool        calculated;   // Normal
   G4bool        validConvex;  // Solid locally convex
   ESide         exitSide;

 public:

   G4ExitNormal(G4ThreeVector norm = G4ThreeVector(0.,0.,0.),
                G4bool        calc = false,
                G4bool        valid= false,
                ESide         side = kNull )
   { exitNormal= norm; calculated= calc; validConvex=valid; exitSide=side;}
};

class G4ReplicaNavigation
{
  public:  // with description

    G4ReplicaNavigation();
   ~G4ReplicaNavigation();

    inline G4bool LevelLocate( G4NavigationHistory& history,
                         const G4VPhysicalVolume *blockedVol,
                         const G4int blockedNum,
                         const G4ThreeVector &globalPoint,
                         const G4ThreeVector* globalDirection,
                         const G4bool pLocatedOnEdge, 
                               G4ThreeVector &localPoint );

    G4double ComputeStep( const G4ThreeVector &globalPoint,
                          const G4ThreeVector &globalDirection,
                          const G4ThreeVector &localPoint,
                          const G4ThreeVector &localDirection,
                          const G4double currentProposedStepLength,
                                G4double &newSafety,
                                G4NavigationHistory &history,
                                G4bool &validExitNormal,
                                G4bool &calculatedExitNormal,
                                G4ThreeVector &exitNormal,
                                G4bool &exiting,
                                G4bool &entering,
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int &blockedReplicaNo );

    G4double ComputeSafety( const G4ThreeVector &globalPoint,
                            const G4ThreeVector &localPoint,
                                  G4NavigationHistory &history,
                            const G4double pProposedMaxLength=DBL_MAX );

    EInside BackLocate( G4NavigationHistory &history,
                  const G4ThreeVector &globalPoint,
                        G4ThreeVector &localPoint,
                  const G4bool &exiting,
                        G4bool &notKnownInside ) const;

    void ComputeTransformation( const G4int replicaNo,
                                      G4VPhysicalVolume *pVol,
                                      G4ThreeVector &point ) const; 
    void ComputeTransformation( const G4int replicaNo,
                                      G4VPhysicalVolume *pVol ) const; 

    EInside Inside( const G4VPhysicalVolume *pVol,
                    const G4int replicaNo,
                    const G4ThreeVector &localPoint ) const;
    G4double DistanceToOut( const G4VPhysicalVolume *pVol,
                            const G4int replicaNo,
                            const G4ThreeVector &localPoint ) const;
    G4double DistanceToOut( const G4VPhysicalVolume *pVol,
                            const G4int replicaNo,
                            const G4ThreeVector &localPoint,
                            const G4ThreeVector &localDirection,
                                  G4ExitNormal& candidateNormal ) const;

    inline G4int GetVerboseLevel() const;
    inline void  SetVerboseLevel(G4int level);
      // Get/Set Verbose(ness) level.
      // [if level>0 && G4VERBOSE, printout can occur]

    inline void  CheckMode(G4bool mode);
      // Run navigation in "check-mode", therefore using additional
      // verifications and more strict correctness conditions.
      // Is effective only with G4VERBOSE set.

  private:

    inline G4int VoxelLocate( const G4SmartVoxelHeader *pHead,
                              const G4ThreeVector &localPoint,
                              const G4int blocked=-1 ) const;

    G4double DistanceToOutPhi( const G4ThreeVector &localPoint,
                               const G4ThreeVector &localDirection,
                               const G4double width,
                               G4ExitNormal& foundNormal ) const;

    G4double DistanceToOutRad( const G4ThreeVector &localPoint,
                               const G4ThreeVector &localDirection,
                               const G4double width,
                               const G4double offset,
                               const G4int replicaNo,
                                     G4ExitNormal& foundNormal ) const;
    inline void SetPhiTransformation( const G4double ang,
                                            G4VPhysicalVolume *pVol=0 ) const;
  private:

    // Invariants - unaltered during navigation
    // **********

    G4bool fCheck; 
    G4int  fVerbose;
      // Configuration parameters

    G4double kCarTolerance, kRadTolerance, kAngTolerance,
             halfkCarTolerance, halfkRadTolerance, halfkAngTolerance,
             fMinStep;
      // Local copy of constants
};

#include "G4ReplicaNavigation.icc"

#endif
