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
//
// $Id: G4PathFinder.hh,v 1.10 2006-05-26 22:00:05 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4Types.hh"
#include <vector>

#include "G4FieldTrack.hh"
// class G4FieldManager;  // #include "G4FieldManager.hh"
// class G4Navigator;
// class G4VPhysicalVolume;
// class G4VCurvedTrajectoryFilter;
class G4TransportationManager; 
class G4Navigator;
#include "G4TouchableHandle.hh"
#include "G4FieldTrack.hh"

enum   ELimited { kDoNot, kUnique, kSharedTransport, kSharedOther, kUndefLimited } ; 

class G4PathFinder
{

 public:  // with description

   static G4PathFinder* GetInstance();  //  Singleton
  ~G4PathFinder();


  //  Attempt next step
  // 
   G4double ComputeStep( const G4FieldTrack      &pFieldTrack,   // Or update non-c
                               G4double           pCurrentProposedStepLength,
			       G4int              navigatorId, 
			       G4int              stepNo,     // See next step / check 
                               G4double          &pNewSafety,   // for this geom 
			       ELimited          &limitedStep, 
			       G4FieldTrack      &EndState );
     // Compute the next geometric Step  -- Curved or linear

   void Locate( const G4ThreeVector& position, 
		const G4ThreeVector& direction,
		G4bool  relativeSearch= true); 
     // Relocate global point in all navigators, and update them.

   void PrepareNewTrack( G4ThreeVector position, G4ThreeVector direction); 
     // Check and cache set of active navigators

   G4TouchableHandle CreateTouchableHandle( G4int navId) const;
   // Also? G4TouchableCreator& GetTouchableCreator( G4int navId ) const; 
   inline G4VPhysicalVolume* GetLocatedVolume( G4int navId ) const; 

  // -----------------------------------------------------------------
   inline void SetChargeMomentumMass( G4double charge,     // in e+ units
                                      G4double momentum,   // in Geant4 units
                                      G4double pMass );  

   inline G4bool         IsParticleLooping() const;

   inline G4int  SetVerboseLevel(G4int lev=-1);
  // inline G4int  GetVerboseLevel() const;

 public:  // with description
   G4double ComputeLinearStep(const G4ThreeVector &pGlobalPoint,
                              const G4ThreeVector &pDirection,
                              G4double pCurrentProposedStepLength,
                              G4double  &pNewSafety,
                              G4bool    &limitedStep, 
                              G4int     stepNo,       // See next step / check 
                              G4int     navId );      // return relevant step
     //  When no field exists or the particle has no charge or EM moment

   inline G4int   GetMaxLoopCount() const;
   inline void    SetMaxLoopCount( G4int new_max );
     // A maximum for the number of steps that a (looping) particle can take.

 public:  // without description
   void MovePoint();  // Signal that location will be moved 

 protected:  // without description
  G4double  DoNextCurvedStep(  const G4FieldTrack  &FieldTrack,
			       G4double            proposedStepLength); 

  G4double  DoNextLinearStep(  const G4FieldTrack  &FieldTrack,
			       G4double            proposedStepLength); 

  void WhichLimited();
  void PrintLimited();   // Print key details out - for debugging
  // 
  // void SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter);
  //   Set the filter that examines & stores 'intermediate' 
  //   curved trajectory points.  Currently only position is stored.

  void ClearState();
  // Clear all the State of this class and its current associates

  inline G4bool UseSafetyForOptimization( G4bool );
      //  Whether to safety to discard 
      //   unneccesary calls to navigator (thus 'optimising' performance)

 protected:
   G4Navigator* GetNavigator(G4int n) const { 
      if( (n>fNoActiveNavigators)||(n<0)){ n=0; }
      return fpNavigator[n]; 
   }

 private:
   G4PathFinder();  //  Singleton 

 private:

  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------
   // std::vector<G4Navigator*> fActiveNavigators; 
     // 
   G4int   fNoActiveNavigators; 
   // G4int   fNoNavigators; 
   G4bool  fNewTrack;               // Flag a new track (ensure first step)

   static const G4int fMaxNav = 8;   // rename to kMaxNoNav ??
   // enum EMaximumNavs { fMaxNav = 8; }

   // Global state (retained during stepping for one track
   G4Navigator*  fpNavigator[fMaxNav];   // G4Navigator** fpNavigator;
   // State after a step computation 
   ELimited      fLimitedStep[fMaxNav];
   G4bool        fLimitTruth[fMaxNav];
   G4double      fCurrentStepSize[fMaxNav]; 
   G4double      fNewSafety[ fMaxNav ];      // Safety for starting point
   G4double      fMinSafety;
   G4double      fMinStep;      // As reported by Navigators -- can be kInfinity
   G4double      fTrueMinStep;  // Corrected in case >= proposed 
   // State after calling 'locate'
   G4VPhysicalVolume* fLocatedVolume[fMaxNav];

   // G4TransportationManager* fpTransportManager; 

   G4FieldTrack    fEndState;
     // End point storage
   G4bool fRelocatedPoint;   //  Signals that point was or is being moved 
                             //  from the position of the last location
                             //   or the endpoint resulting from ComputeStep 
                             //   -- invalidates fEndState

   G4int           fLastStepNo, fCurrentStepNo; 
   G4bool      fParticleIsLooping;
   G4int  fVerboseLevel;
     // For debuging purposes
   G4int  fMax_loop_count;
    // Limit for the number of sub-steps taken in one call to ComputeStep
};

// ********************************************************************
// Inline methods.
// ********************************************************************
inline G4VPhysicalVolume* G4PathFinder::GetLocatedVolume( G4int navId ) const
{  G4VPhysicalVolume* vol=0;  
   if( (navId < fMaxNav) && (navId >=0) ) { vol= fLocatedVolume[navId]; }
   return vol; 
}

inline G4int  G4PathFinder::SetVerboseLevel(G4int newLevel)
{  G4int old= fVerboseLevel;  fVerboseLevel= newLevel; return old;
}
// #include "G4PathFinder.icc"

#endif 
