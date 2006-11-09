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
// $Id: G4MultiNavigator.hh,v 1.1 2006-11-09 13:55:29 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MultiNavigator
//
// Class description:
//
// Polls the navigators of several geometries to identify the next 
// boundary. 
//
// History:
// - Created.                                John Apostolakis  Nov 2006
// *********************************************************************

#ifndef G4MULTINAVIGATOR_HH
#define G4MULTINAVIGATOR_HH

#include "geomdefs.hh"
#include "G4ThreeVector.hh"
// #include "G4Types.hh"

#include "G4Navigator.hh"

// #include "G4AffineTransform.hh"
// #include "G4RotationMatrix.hh"

#include "G4TouchableHistoryHandle.hh"

#include "G4NavigationHistory.hh"

#include "G4Elimited.hh"
class G4TransportationManager; 

#include <iostream>

class G4VPhysicalVolume;

class G4MultiNavigator : public G4Navigator
{
  public:  // with description

  friend std::ostream& operator << (std::ostream &os, const G4Navigator &n);

  G4MultiNavigator();
    // Constructor - initialisers and setup.

  ~G4MultiNavigator();
    // Destructor. No actions.

  // virtual 
  G4double ComputeStep(const G4ThreeVector &pGlobalPoint,
		       const G4ThreeVector &pDirection,
		       const G4double      pCurrentProposedStepLength,
		             G4double      &pNewSafety);
    // Returns the distance to the next boundary of any geometry

  G4double ObtainFinalStep( G4int        navigatorId, 
			    G4double     &pNewSafety,     // for this geom 
			    G4double     &minStep,
			    ELimited     &limitedStep); 
    // Get values for a single geometry

  void PrepareNavigators(); 
    //  Find which geometries are registered for this particles - and keep info
  void PrepareNewTrack( const G4ThreeVector position, 
			const G4ThreeVector direction ); 
    //  Prepare Navigators and locate 

  // virtual
  G4VPhysicalVolume* ResetHierarchyAndLocate(const G4ThreeVector &point,
                                             const G4ThreeVector &direction,
                                             const G4TouchableHistory &h);

    // Resets the geometrical hierarchy for all geometries.
    //  Use the touchable history for the first (mass) geometry.
    //  Return the volume in the first (mass) geometry.
    // 
    // Important Note: In order to call this the geometries MUST be closed.

  // virtual
  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                             const G4ThreeVector* direction=0,
                                             const G4bool pRelativeSearch=true,
                                             const G4bool ignoreDirection=true);
    // Locate in all geometries.
    // Return the volume in the first (mass) geometry
    //  Maintain vector of other volumes,  to be returned separately
    // 
    // Important Note: In order to call this the geometry MUST be closed.

  // virtual
  void LocateGlobalPointWithinVolume(const G4ThreeVector& position);
    // Relocate in all geometries for point that has not changed volume
    // (ie is within safety  in all geometries or is distance less that 
    // along the direction of a computed step.

  // void SetGeometricallyLimitedStep();
    // Inform the navigator that the previous Step calculated
    // by the geometry was taken in its entirety.

  // virtual 
  G4double ComputeSafety(const G4ThreeVector &globalpoint,
                                 const G4double pProposedMaxLength = DBL_MAX);
    // Calculate the isotropic distance to the nearest boundary 
    // in any geometry from the specified point in the global coordinate system. 
        // The geometry must be closed.

  // virtual 
  G4TouchableHistoryHandle CreateTouchableHistoryHandle() const;
    // Returns a reference counted handle to a touchable history.

  // virtual 
  //   G4ThreeVector GetLocalExitNormal(G4bool* valid);
  //  ---> Not overloading currently   TODO: check if relevant

  // Check relevance, use
  // G4bool EnteredDaughterVolume() const;
  // G4bool ExitedMotherVolume() const;
  // void  CheckMode(G4bool mode);

  // void ResetStackAndState();
  // G4int SeverityOfZeroStepping( G4int* noZeroSteps ) const; 

  G4Navigator* GetNavigator(G4int n) const { 
     if( (n>fNoActiveNavigators)||(n<0)){ n=0; }
     return fpNavigator[n]; 
  }

 protected:  // with description

  // virtual 
  void ResetState();
    // Utility method to reset the navigator state machine.

  // virtual 
  void SetupHierarchy();
    // Renavigate & reset hierarchy described by current history
    // o Reset volumes
    // o Recompute transforms and/or solids of replicated/parameterised
    //   volumes.

  void WhichLimited(); // Flag which processes limited the step
  void PrintLimited(); // Auxiliary, debugging printing 
  void CheckMassWorld(); 

 private:

   // STATE Information 

   G4int   fNoActiveNavigators; 
   static const G4int fMaxNav = 8;   // rename to kMaxNoNav ??
   G4VPhysicalVolume* fLastMassWorld; 

   // Global state (retained during stepping for one track
   G4Navigator*  fpNavigator[fMaxNav];   // G4Navigator** fpNavigator;

   // State after a step computation 
   ELimited      fLimitedStep[fMaxNav];
   G4bool        fLimitTruth[fMaxNav];
   G4double      fCurrentStepSize[fMaxNav]; 
   G4double      fNewSafety[ fMaxNav ];      // Safety for starting point

   // Lowest values - determine step length, and safety 
   G4double      fMinStep;      // As reported by Navigators -- can be kInfinity
   G4double      fMinSafety;
   G4double      fTrueMinStep;  // Corrected in case fMinStep >= proposed 

   // State after calling 'locate'
   G4VPhysicalVolume* fLocatedVolume[fMaxNav];
   G4ThreeVector      fLastLocatedPosition; 

   // cache of safety information 
   G4ThreeVector fSafetyLocation;       //  point where ComputeSafety is called
   G4double      fMinSafety_atSafLocation; // /\ corresponding value of safety
   G4ThreeVector fPreStepLocation;      //  point where last ComputeStep called
   G4double      fMinSafety_PreStepPt;  //   /\ corresponding value of safety

   // STATE Auxiliary variables 
   G4int  fVerboseLevel;
   G4TransportationManager* pTransportManager; // Cache for frequent use
};

// #include "G4MultiNavigator.icc"

#endif

