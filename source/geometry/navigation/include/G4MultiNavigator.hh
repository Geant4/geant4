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
// $Id: G4MultiNavigator.hh 87869 2015-01-16 08:24:36Z gcosmo $
//
//
// class G4MultiNavigator
//
// Class description:
//
// Utility class for polling the navigators of several geometries to
// identify the next  boundary. 

// History:
// - Created. John Apostolakis, November 2006
// *********************************************************************

#ifndef G4MULTINAVIGATOR_HH
#define G4MULTINAVIGATOR_HH

#include <iostream>

#include "geomdefs.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"

#include "G4TouchableHistoryHandle.hh"

#include "G4NavigationHistory.hh"

enum  ELimited { kDoNot,kUnique,kSharedTransport,kSharedOther,kUndefLimited };

class G4TransportationManager;
class G4VPhysicalVolume;

class G4MultiNavigator : public G4Navigator
{
  public:  // with description

  friend std::ostream& operator << (std::ostream &os, const G4Navigator &n);

  G4MultiNavigator();
    // Constructor - initialisers and setup.

  ~G4MultiNavigator();
    // Destructor. No actions.

  G4double ComputeStep(const G4ThreeVector &pGlobalPoint,
                       const G4ThreeVector &pDirection,
                       const G4double      pCurrentProposedStepLength,
                             G4double      &pNewSafety);
    // Return the distance to the next boundary of any geometry

  G4double ObtainFinalStep( G4int        navigatorId, 
                            G4double     &pNewSafety,     // for this geom 
                            G4double     &minStepLast,
                            ELimited     &limitedStep); 
    // Get values for a single geometry

  void PrepareNavigators(); 
    // Find which geometries are registered for this particles, and keep info
  void PrepareNewTrack( const G4ThreeVector position, 
                        const G4ThreeVector direction ); 
    // Prepare Navigators and locate 

  G4VPhysicalVolume* ResetHierarchyAndLocate(const G4ThreeVector &point,
                                             const G4ThreeVector &direction,
                                             const G4TouchableHistory &h);
    // Reset the geometrical hierarchy for all geometries.
    // Use the touchable history for the first (mass) geometry.
    // Return the volume in the first (mass) geometry.
    // 
    // Important Note: In order to call this the geometries MUST be closed.

  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                         const G4ThreeVector* direction=0,
                                         const G4bool pRelativeSearch=true,
                                         const G4bool ignoreDirection=true);
    // Locate in all geometries.
    // Return the volume in the first (mass) geometry
    //  Maintain vector of other volumes,  to be returned separately
    // 
    // Important Note: In order to call this the geometry MUST be closed.

   void LocateGlobalPointWithinVolume(const G4ThreeVector& position); 
    // Relocate in all geometries for point that has not changed volume
    // (ie is within safety  in all geometries or is distance less that 
    // along the direction of a computed step.

  G4double ComputeSafety(const G4ThreeVector &globalpoint,
                         const G4double pProposedMaxLength = DBL_MAX,
                         const G4bool keepState = false);
    // Calculate the isotropic distance to the nearest boundary 
    // in any geometry from the specified point in the global coordinate
    // system. The geometry must be closed.

  G4TouchableHistoryHandle CreateTouchableHistoryHandle() const;
    // Returns a reference counted handle to a touchable history.

  virtual G4ThreeVector GetLocalExitNormal(G4bool* obtained); // const
  virtual G4ThreeVector GetLocalExitNormalAndCheck(const G4ThreeVector &CurrentE_Point,
						   G4bool* obtained); // const
  virtual G4ThreeVector GetGlobalExitNormal(const G4ThreeVector &CurrentE_Point,
 					           G4bool* obtained); // const
    // Return Exit Surface Normal and validity too.
    // Can only be called if the Navigator's last Step either
    //  - has just crossed a volume geometrical boundary and relocated, or
    //  - has arrived at a boundary in a ComputeStep
    // It returns the Normal to the surface pointing out of the volume that
    //   was left behind and/or into the volume that was entered.
    // Convention:x
    //   The *local* normal is in the coordinate system of the *final* volume.
    // Restriction:
    //   Normals are not available for replica volumes (returns obtained= false)

 public:  // without description

  G4Navigator* GetNavigator(G4int n) const
  { 
     if( (n>fNoActiveNavigators)||(n<0)){ n=0; }
     return fpNavigator[n]; 
  }

 protected:  // with description

  void ResetState();
    // Utility method to reset the navigator state machine.

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
   static const G4int fMaxNav = 16;   // rename to kMaxNoNav ??
   G4VPhysicalVolume* fLastMassWorld; 

   // Global state (retained during stepping for one track
   G4Navigator*  fpNavigator[fMaxNav];   // G4Navigator** fpNavigator;

   // State after a step computation 
   ELimited      fLimitedStep[fMaxNav];
   G4bool        fLimitTruth[fMaxNav];
   G4double      fCurrentStepSize[fMaxNav]; 
   G4double      fNewSafety[ fMaxNav ];      // Safety for starting point
   G4int         fNoLimitingStep;       // How many geometries limited the step
   G4int         fIdNavLimiting;        // Id of Navigator limiting step (if only one limits)

   // Lowest values - determine step length, and safety 
   G4double      fMinStep;      // As reported by Navigators. Can be kInfinity
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

   G4TransportationManager* pTransportManager; // Cache for frequent use
};

#endif
