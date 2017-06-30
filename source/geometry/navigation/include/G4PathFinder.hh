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
// $Id: G4PathFinder.hh 103219 2017-03-22 11:30:15Z gcosmo $
// 
// class G4PathFinder 
//
// Class description:
// 
// This class directs the lock-stepped propagation of a track in the 
// 'mass' and other parallel geometries.  It ensures that tracking 
// in a magnetic field sees these parallel geometries at each trial step, 
// and that the earliest boundary limits the step.
// 
// For the movement in field, it relies on the class G4PropagatorInField
//
// History:
// -------
//  7.10.05 John Apostolakis,  Draft design 
// 26.04.06 John Apostolakis,  Revised design and first implementation 
// ---------------------------------------------------------------------------
#ifndef G4PATHFINDER_HH 
#define G4PATHFINDER_HH  1

#include <vector>
#include "G4Types.hh"

#include "G4FieldTrack.hh"

class G4TransportationManager; 
class G4Navigator;

#include "G4TouchableHandle.hh"
#include "G4FieldTrack.hh"
#include "G4MultiNavigator.hh"

class G4PropagatorInField;

class G4PathFinder
{

 public:  // with description

   static G4PathFinder* GetInstance();
     //
     // Retrieve singleton instance and create it if not existing.

   static G4PathFinder* GetInstanceIfExist();
     //
     // Retrieve singleton instance pointer.

   G4double ComputeStep( const G4FieldTrack &pFieldTrack,
                         G4double  pCurrentProposedStepLength,
                         G4int     navigatorId, // Identifies the geometry
                         G4int     stepNo,      // See next step/check 
                         G4double &pNewSafety,  // Only for this geometry
                         ELimited &limitedStep,      
                         G4FieldTrack       &EndState, 
                         G4VPhysicalVolume*  currentVolume );
     //
     // Compute the next geometric Step  -- Curved or linear
     //   If it is called with a larger 'stepNo' it will execute a new step;
     //   if 'stepNo' is same as last call, then the results for 
     //   the geometry with Id. number 'navigatorId' will be returned. 

   void Locate( const G4ThreeVector& position, 
                const G4ThreeVector& direction,
                      G4bool  relativeSearch=true); 
     //
     // Make primary relocation of global point in all navigators,
     // and update them.

   void ReLocate( const G4ThreeVector& position ); 
     //
     // Make secondary relocation of global point (within safety only) 
     // in all navigators, and update them.

   void PrepareNewTrack( const G4ThreeVector& position,
                         const G4ThreeVector& direction,
                               G4VPhysicalVolume* massStartVol=0); 
     //
     // Check and cache set of active navigators.

   void EndTrack();
     // Signal end of tracking of current track.  
     //   Reset internal state
     //   Inform TransportationManager to use 'ordinary' Navigator

    G4TouchableHandle CreateTouchableHandle( G4int navId ) const;
    inline G4VPhysicalVolume* GetLocatedVolume( G4int navId ) const;

    G4bool RecheckDistanceToCurrentBoundary(
                                          const G4ThreeVector &pGlobalPoint,
                                          const G4ThreeVector &pDirection,
                                          const G4double pCurrentProposedStepLength,
                                          G4double  *prDistance,
                                          G4double  *prNewSafety= 0)const;
   // Trial method for checking potential displacement for MS

   // -----------------------------------------------------------------
  
   inline G4bool   IsParticleLooping() const;

   inline G4double GetCurrentSafety() const;
     // Minimum value of safety after last ComputeStep
   inline G4double GetMinimumStep() const;      
     // Get the minimum step size from the last ComputeStep call
     //   - in case full step is taken, this is kInfinity
   inline unsigned int  GetNumberGeometriesLimitingStep() const; 

   G4double ComputeSafety( const G4ThreeVector& globalPoint); 
     // Recompute safety for the relevant point the endpoint of the last step!!
     // Maintain vector of individual safety values (for next method)

   G4double ObtainSafety( G4int navId, G4ThreeVector& globalCenterPoint );
     // Obtain safety for navigator/geometry navId for last point 'computed'
     //   --> last point for which ComputeSafety was called
     //   Returns the point (center) for which this safety is valid

   void EnableParallelNavigation( G4bool enableChoice=true ); 
     //
     // Must call it to ensure that PathFinder is prepared,  
     // especially for curved tracks. If true it switches PropagatorInField
     // to use MultiNavigator. Must call it with false to undo (=PiF use
     // Navigator for tracking!)

   inline G4int  SetVerboseLevel(G4int lev=-1);

 public:  // with description

   inline G4int   GetMaxLoopCount() const;
   inline void    SetMaxLoopCount( G4int new_max );
     //
     // A maximum for the number of steps that a (looping) particle can take.

 public:  // without description

   inline void MovePoint();
     //
     // Signal that location will be moved -- internal use primarily

   // To provide best compatibility between Coupled and Old Transportation
   //   the next two methods are provided:
   G4double LastPreSafety( G4int navId, G4ThreeVector& globalCenterPoint, G4double& minSafety ); 
     // Obtain last safety needed in ComputeStep (for geometry navId)
     //   --> last point at which ComputeStep recalculated safety
     //   Returns the point (center) for which this safety is valid
     //    and also the minimum safety over all navigators (ie full)

   void PushPostSafetyToPreSafety(); 
     // Tell PathFinder to copy PostStep Safety to PreSafety (for use at next step)

   G4String& LimitedString( ELimited lim );
     // Convert ELimited to string

   ~G4PathFinder();
     // Destructor

 protected:  // without description

  G4double  DoNextLinearStep(  const G4FieldTrack  &FieldTrack,
                               G4double            proposedStepLength); 

  G4double  DoNextCurvedStep(  const G4FieldTrack  &FieldTrack,
                               G4double            proposedStepLength,
                               G4VPhysicalVolume*  pCurrentPhysVolume); 

  void WhichLimited();
  void PrintLimited();
  //
  // Print key details out - for debugging

  // void ClearState();
  //
  // Clear all the State of this class and its current associates

  inline G4bool UseSafetyForOptimization( G4bool );
  //
  // Whether use safety to discard unneccesary calls to navigator

  void ReportMove( const G4ThreeVector& OldV, const G4ThreeVector& NewV, const G4String& Quantity ) const; 
  // Helper method to report movement (likely of initial point)

 protected:

   G4PathFinder();  //  Singleton 

  inline G4Navigator* GetNavigator(G4int n) const;

 private:

  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------

   G4MultiNavigator *fpMultiNavigator; 
   //
   //  Object that enables G4PropagatorInField to see many geometries

   G4int   fNoActiveNavigators; 
   G4bool  fNewTrack;               // Flag a new track (ensure first step)

   static const G4int fMaxNav = 16;  // rename to kMaxNoNav ??

   // Global state (retained during stepping for one track)

   G4Navigator*  fpNavigator[fMaxNav];

   // State changed in a step computation

   ELimited      fLimitedStep[fMaxNav];
   G4bool        fLimitTruth[fMaxNav];
   G4double      fCurrentStepSize[fMaxNav]; 
   G4int         fNoGeometriesLimiting;  //  How many processes contribute to limit

   G4ThreeVector fPreSafetyLocation;    //  last initial position for which safety evaluated
   G4double      fPreSafetyMinValue;    //   /\ corresponding value of full safety
   G4double      fPreSafetyValues[ fMaxNav ]; //   Safeties for the above point
   // This part of the state can be retained for severall calls --> CARE

   G4ThreeVector fPreStepLocation;      //  point where last ComputeStep called
   G4double      fMinSafety_PreStepPt;  //   /\ corresponding value of full safety
   G4double      fCurrentPreStepSafety[ fMaxNav ]; //   Safeties for the above point
   // This changes at each step, 
   //   so it can differ when steps inside min-safety are made

   G4bool        fPreStepCenterRenewed;   // Whether PreSafety coincides with PreStep point 

   G4double      fMinStep;      // As reported by Navigators -- can be kInfinity
   G4double      fTrueMinStep;  // Corrected in case >= proposed

   // State after calling 'locate'

   G4VPhysicalVolume* fLocatedVolume[fMaxNav];
   G4ThreeVector      fLastLocatedPosition; 

   // State after calling 'ComputeStep' (others member variables will be affected)
   G4FieldTrack    fEndState;           // Point, velocity, ... at proposed step end
   G4bool          fFieldExertedForce;  // In current proposed step

   G4bool fRelocatedPoint;   //  Signals that point was or is being moved 
                             //  from the position of the last location
                             //   or the endpoint resulting from ComputeStep 
                             //   -- invalidates fEndState

   // State for 'ComputeSafety' and related methods
   G4ThreeVector fSafetyLocation;       //  point where ComputeSafety is called
   G4double      fMinSafety_atSafLocation; // /\ corresponding value of safety
   G4double      fNewSafetyComputed[ fMaxNav ];  // Safeties for last ComputeSafety

   // State for Step numbers 
   G4int         fLastStepNo, fCurrentStepNo; 

   G4int         fVerboseLevel;            // For debuging purposes

   G4TransportationManager* fpTransportManager; // Cache for frequent use
   G4PropagatorInField* fpFieldPropagator;

   G4double kCarTolerance;

   static G4ThreadLocal G4PathFinder* fpPathFinder;
};

// ********************************************************************
// Inline methods.
// ********************************************************************

inline G4VPhysicalVolume* G4PathFinder::GetLocatedVolume( G4int navId ) const
{
  G4VPhysicalVolume* vol=0;  
  if( (navId < fMaxNav) && (navId >=0) ) { vol= fLocatedVolume[navId]; }
  return vol; 
}

inline G4int  G4PathFinder::SetVerboseLevel(G4int newLevel)
{
  G4int old= fVerboseLevel;  fVerboseLevel= newLevel; return old;
}

inline G4double G4PathFinder::GetMinimumStep() const
{ 
  return fMinStep; 
} 

inline unsigned int G4PathFinder::GetNumberGeometriesLimitingStep() const
{
  unsigned int noGeometries=fNoGeometriesLimiting;
  return noGeometries; 
}

inline G4double G4PathFinder::GetCurrentSafety() const
{
  return fMinSafety_PreStepPt;
}

inline void G4PathFinder::MovePoint()
{
  fRelocatedPoint= true;
}

inline G4Navigator* G4PathFinder::GetNavigator(G4int n) const
{ 
  if( (n>fNoActiveNavigators)||(n<0)) { n=0; }
  return fpNavigator[n];
}

inline G4double      G4PathFinder::ObtainSafety( G4int navId, G4ThreeVector& globalCenterPoint )
{
  globalCenterPoint= fSafetyLocation; 
  //  navId = std::min( navId, fMaxNav-1 ); 
  return  fNewSafetyComputed[ navId ];
}

inline G4double  G4PathFinder::LastPreSafety( G4int navId, 
					      G4ThreeVector& globalCenterPoint, 
					      G4double& minSafety )
{
  globalCenterPoint= fPreSafetyLocation;
  minSafety=         fPreSafetyMinValue;
  //  navId = std::min( navId, fMaxNav-1 ); 
  return  fPreSafetyValues[ navId ];
}
#endif 
