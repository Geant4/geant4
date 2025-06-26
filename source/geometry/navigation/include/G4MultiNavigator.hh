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
// G4MultiNavigator
//
// Class description:
//
// Utility class for polling the navigators of several geometries to
// identify the next boundary. 

// Author: John Apostolakis (CERN), November 2006
// --------------------------------------------------------------------
#ifndef G4MULTINAVIGATOR_HH
#define G4MULTINAVIGATOR_HH 1

#include <iostream>

#include "geomdefs.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"

#include "G4TouchableHandle.hh"

#include "G4NavigationHistory.hh"

enum  ELimited { kDoNot,kUnique,kSharedTransport,kSharedOther,kUndefLimited };

class G4TransportationManager;
class G4VPhysicalVolume;

/**
 * @brief G4MultiNavigator is a utility class for polling the navigators
 * of several geometries to identify the next boundary.
 */

class G4MultiNavigator : public G4Navigator
{
  public:

    friend std::ostream& operator << (std::ostream& os, const G4Navigator& n);

    /**
     * Constructor and default Destructor.
     */
    G4MultiNavigator();
    ~G4MultiNavigator() override = default;

    /**
     * Computes the distance to the next boundary of any geometry.
     *  @param[in] pGlobalPoint The point in global coordinates system.
     *  @param[in] pDirection The normalised vector direction.
     *  @param[in] pCurrentProposedStepLength Current proposed step length.
     *  @param[in,out] newSafety New safety.
     *  @returns Length from current point to next boundary surface along
     *           @p pDirection.
     */
    G4double ComputeStep( const G4ThreeVector& pGlobalPoint,
                          const G4ThreeVector& pDirection,
                          const G4double       pCurrentProposedStepLength,
                                G4double&      pNewSafety ) override;

    /**
     * Gets values for a single geometry.
     *  @param[in] navigatorId The navigator identifier.
     *  @param[in,out] pnewSafety New safety for this geometry.
     *  @param[in,out] minStepLast The last minimum step returned.
     *  @param[in,out] limitedStep The step characterisation returned.
     *  @returns The step size for the geometry associated to 'navigatorId'.
     */
    G4double ObtainFinalStep( G4int        navigatorId, 
                              G4double&    pNewSafety,     // for this geom 
                              G4double&    minStepLast,
                              ELimited&    limitedStep ); 

    /**
     * Finds which geometries are registered for this particles, and keeps info.
     */
    void PrepareNavigators(); 

    /**
     * Prepares Navigators and locates.
     *  @param[in] position The position point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     */
    void PrepareNewTrack( const G4ThreeVector& position, 
                          const G4ThreeVector direction ); 

    /**
     * Resets the geometrical hierarchy for all geometries.
     * Use the touchable history for the first (mass) geometry.
     *  @note In order to call this the geometries MUST be closed.
     *  @param[in] point The  point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     *  @param[in] h The touchable history to be used for initialisation.
     *  @returns The pointer to the volume in the first (mass) geometry.
     */
    G4VPhysicalVolume* ResetHierarchyAndLocate( const G4ThreeVector& point,
                                     const G4ThreeVector& direction,
                                     const G4TouchableHistory& h ) override;

    /**
     * Locates the point in all geometries.
     * Maintains a vector of other volumes, to be returned separately.
     *  @note In order to call this the geometry MUST be closed.
     *  @param[in] point The point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     *  @param[in] pRelativeSearch Flag to specify where search starts from.
     *  @param[in] ignoreDirection Flag to specify if to use direction or not.
     *  @returns The volume in the first (mass) geometry.
     */
    G4VPhysicalVolume* LocateGlobalPointAndSetup( const G4ThreeVector& point,
                                      const G4ThreeVector* direction = nullptr,
                                      const G4bool pRelativeSearch = true,
                                      const G4bool ignoreDirection = true) override;

    /**
     * Relocates in all geometries for point that has not changed volume,
     * i.e. is within safety in all geometries or its distance is less that 
     * along the direction of a computed step.
     *  @param[in] position The position point in global coordinates system.
     */
    void LocateGlobalPointWithinVolume( const G4ThreeVector& position ) override; 

    /**
     * Calculates the isotropic distance to the nearest boundary in any
     * geometry from the specified point in the global coordinates system.
     *  @note The geometry must be closed.
     *  @param[in] globalpoint The point in global coordinates system.
     *             The point must be within the current volume.
     *  @param[in] pProposedMaxLength The proposed maximum length is used
     *             to avoid volume safety calculations.
     *  @param[in] keepState Flag to instruct keeping the state (default false)
     *             to ensure minimum side effects from the call.
     *  @returns Length from current point to closest boundary surface.
     *           The value returned is usually an underestimate.  
     */
    G4double ComputeSafety( const G4ThreeVector& globalpoint,
                            const G4double pProposedMaxLength = DBL_MAX,
                            const G4bool keepState = false ) override;

    /**
     * Returns a reference counted handle to a touchable history.
     */
    G4TouchableHandle CreateTouchableHistoryHandle() const override;

    /**
     * Obtains the Normal vector to a surface (in local coordinates)
     * pointing out of previous volume and into current volume
     * Convention: the *local* normal is in the coordinate system of the
     * *final* volume. The method takes full care about how to calculate
     * this normal, but if the surfaces are not convex it will return
     * valid=false.
     *  @param[in,out] obtained Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    G4ThreeVector GetLocalExitNormal( G4bool* obtained ) override;

    /**
     * Obtains the Normal vector to a surface (in local coordinates)
     * pointing out of previous volume and into current volume, and
     * checks the current point against expected 'local' value.
     * Convention: the *local* normal is in the coordinate system of the
     * *final* volume. The method takes full care about how to calculate
     * this normal, but if the surfaces are not convex it will return
     * valid=false.
     *  @param[in] point Point in global coordinates system to compare to.
     *  @param[in,out] obtained Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    G4ThreeVector GetLocalExitNormalAndCheck( const G4ThreeVector& point,
                                              G4bool* obtained ) override;

    /**
     * Obtains the Normal vector to a surface (in global coordinates)
     * pointing out of previous volume and into current volume
     * The method takes full care about how to calculate the normal,
     * but if the surfaces are not convex it will return valid=false.
     *  @param[in] point Point in global coordinates system to compare to.
     *  @param[in,out] obtained Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    G4ThreeVector GetGlobalExitNormal( const G4ThreeVector& point,
                                             G4bool* obtained ) override;

    /**
     * Returns a pointer to a navigator, given its index.
     */
    inline G4Navigator* GetNavigator( G4int n ) const;

  protected:

    /**
     * Utility method to reset the navigator state machine.
     */
    void ResetState() override;

    /**
     * Renavigates & resets hierarchy described by the current history,
     * i.e. resets volumes and recomputes transforms and/or solids of
     * replicated/parameterised volumes.
     */
    void SetupHierarchy() override;

    /**
     * Flags which processes limited the step.
     */
    void WhichLimited(); 

    /**
     * Auxiliary, debugging printing.
     */
    void PrintLimited();

    /**
     * Checks if mass world pointed has been changed => issues and exception.
     */
    void CheckMassWorld(); 

  private:

    // STATE Information 

    G4int fNoActiveNavigators = 0; 
    static const G4int fMaxNav = 16;
    G4VPhysicalVolume* fLastMassWorld = nullptr; 

    /** Global state (retained during stepping for one track). */
    G4Navigator* fpNavigator[fMaxNav];

    // State after a step computation 
    //
    ELimited fLimitedStep[fMaxNav];
    G4bool   fLimitTruth[fMaxNav];
    G4double fCurrentStepSize[fMaxNav]; 
    G4double fNewSafety[ fMaxNav ]; // Safety for starting point
    G4int    fNoLimitingStep = -1;  // How many geometries limited the step
    G4int    fIdNavLimiting = -1;   // Id of Navigator limiting step

    // Lowest values - determine step length, and safety
    //
    G4double fMinStep = -kInfinity;     // As reported by Navigators
    G4double fMinSafety = -kInfinity;
    G4double fTrueMinStep = -kInfinity; // Corrected if fMinStep>=proposed

    // State after calling 'locate'
    //
    G4VPhysicalVolume* fLocatedVolume[fMaxNav];
    G4ThreeVector      fLastLocatedPosition; 

    // Cache of safety information
    //
    G4ThreeVector fSafetyLocation; // point where ComputeSafety() is called
    G4double fMinSafety_atSafLocation = -1.0; // - corresponding value of safety
    G4ThreeVector fPreStepLocation; // point where last ComputeStep() called
    G4double fMinSafety_PreStepPt = -1.0; // - corresponding value of safety

    G4TransportationManager* pTransportManager; // Cache for frequent use
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4Navigator* G4MultiNavigator::GetNavigator( G4int n ) const
{ 
  if( (n>fNoActiveNavigators) || (n<0) ) { n=0; }
  return fpNavigator[n]; 
}

#endif
