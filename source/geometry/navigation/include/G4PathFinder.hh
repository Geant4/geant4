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
// $Id: G4PathFinder.hh,v 1.1 2006-04-26 13:21:43 japost Exp $
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
// Key Method:
//              ComputeStep(..)
// History:
// -------
//  7.10.05 John Apostolakis,  design and implementation 
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

class G4PathFinder
{

 public:  // with description

   GetInstance();  //  Singleton
  ~G4PathFinder();

  //  Attempt next step
  // 
   G4double ComputeStep( const G4FieldTrack      &pFieldTrack,   // Or update non-c
                               G4double           pCurrentProposedStepLength,
			       G4int              navigatorId, 
			       G4int              stepNo,     // See next step / check 
                               G4double          &pNewSafety,   // for this geom 
			       G4bool            &limitedStep, 
			       G4FieldTrack      &EndState );
     // Compute the next geometric Step  -- Curved or linear

   void PrepareNewTrack( G4ThreeVector position, G4ThreeVector direction); 
     // Check and cache set of active navigators

   G4TouchableHistoryHandle CreateTouchableHandle( navId ) const;
   // Also? G4TouchableCreator& GetTouchableCreator( navId ) const; 

  // -----------------------------------------------------------------
   inline void SetChargeMomentumMass( G4double charge,     // in e+ units
                                      G4double momentum,   // in Geant4 units
                                      G4double pMass );  

   inline G4bool         IsParticleLooping() const;

   inline G4int  SetVerboseLevel( G4int verbose= -1 );
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

 protected:  // without description
  // 
  void SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter);
  // Set the filter that examines & stores 'intermediate' 
  //  curved trajectory points.  Currently only position is stored.

  void ClearState();
  // Clear all the State of this class and its current associates

  inline G4bool UseSafetyForOptimization( G4bool );
      //  Whether to safety to discard 
      //   unneccesary calls to navigator (thus 'optimising' performance)

 protected:  // with description

 private:
   G4PathFinder();  //  Singleton 

 private:

  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------
   std::vec<G4Navigator> fActiveNavigators; 
     // 
   G4int   fNoActiveNavigators; 
   G4int   fNoTotalNavigators; 

    
   G4double fLimitedStep[MaxNav]; 
   G4double fCurrentStepSize[MaxNav]; 

   G4FieldTrack    End_PointAndTangent;
     // End point storage
   G4bool      fParticleIsLooping;
   G4int  fVerboseLevel;
     // For debuging purposes
   G4int  fMax_loop_count;
    // Limit for the number of sub-steps taken in one call to ComputeStep


};

// ********************************************************************
// Inline methods.
// ********************************************************************
#include "G4PathFinder.icc"

#endif 


#if 0

// Methods moved to G4NavigatorDispatcher (or similar class)

  // Initialisation, activation, de-activation methods
  // --------------------------
   G4int   RegisterNavigator( G4Navigator *pNavigator ); 
   G4int   DeregisterNavigator( G4Navigator *pNavigator ); 
     // Registers the navigator and returns its ID
     //   -> registering a new navigator will clear all Active flags 
     //       so it can only be done before tracking starts.
   G4bool  ActivateNavigator( G4int navId );
   G4bool  DeactivateNavigator( G4int navId ); 
     // Activate (de-activate) the relevant navigator (return true=ok) 
     //  First one to be activated will be 'lead' navigator (its call will initiative step)
   ClearActiveNavigators(); 
     // Clears active flag for all navigators - must be called at the end of a track

#endif


#if 0
  // Proposal 2 or internal method ??
  G4double  DoNextCurvedStep(  G4FieldTrack  &pFieldTrack,
			       G4double       pCurrentProposedStepLength,
                               G4double      &pOverallNewSafety, 
			       G4double       stepNo ); 
     // Initiates the curved step for all geometries

  G4double  CheckCurvedStep(   G4int         navId 
			       G4double      &pNewSafety, 
                               G4bool        &limitedStep );
     // Retrieves the step information for navigator 'navId'

endif
