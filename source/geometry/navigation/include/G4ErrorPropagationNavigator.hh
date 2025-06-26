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
// G4ErrorPropagationNavigator
//
// Class Description:
//
// Class for performing double navigation in the detector geometry and
// on the target surface for error propagation. It overloads ComputeStep()
// and ComputeSafety() methods.

// Author: Pedro Arce (CIEMAT), September 2004
// --------------------------------------------------------------------
#ifndef G4ErrorPropagationNavigator_hh
#define G4ErrorPropagationNavigator_hh 1

#include "G4Navigator.hh"
#include "G4ThreeVector.hh"

/**
 * @brief G4ErrorPropagationNavigator is a class for performing double
 * navigation in the detector geometry and on the target surface for error
 * propagation. It overloads ComputeStep() and ComputeSafety() methods.
 */

class G4ErrorPropagationNavigator : public G4Navigator
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4ErrorPropagationNavigator() = default;
    ~G4ErrorPropagationNavigator() override = default;
  
    /**
     * Calls the navigation in the detector geometry and then checks
     * if the distance to surface is smaller than the proposed step.
     *  @param[in] pGlobalPoint The point in global coordinates system.
     *  @param[in] pDirection The normalised vector direction.
     *  @param[in] pCurrentProposedStepLength Current proposed step length.
     *  @param[in,out] newSafety New safety.
     *  @returns Length from current point to next boundary surface along
     *           @p pDirection.
     */
    G4double ComputeStep (const G4ThreeVector& pGlobalPoint,
                          const G4ThreeVector& pDirection,
                          const G4double pCurrentProposedStepLength,
                                G4double &pNewSafety) override;

    /**
     * Calls the navigation in the detector geometry and then checks
     * if the distance to surface is smaller than the proposed safety.
     *  @param[in] globalpoint The point in global coordinates system.
     *             The point must be within the current volume.
     *  @param[in] pProposedMaxLength The proposed maximum length is used
     *             to avoid volume safety calculations.
     *  @param[in] keepState Flag to instruct keeping the state (default true)
     *             to ensure minimum side effects from the call.
     *  @returns Length from current point to closest boundary surface.
     *           The value returned is usually an underestimate.  
     */
    G4double ComputeSafety(const G4ThreeVector& globalpoint,
                           const G4double pProposedMaxLength = DBL_MAX,
                           const G4bool keepState = true) override;
  
    /**
     * Returns Exit Surface Normal and validity too. Can only be called if
     * the Navigator's last Step has crossed a volume geometrical boundary.
     * Normal points out of the volume exited and/or into the volume entered.
     *  @param[in] point Point in global coordinates system to compare to.
     *  @param[in,out] valid Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    G4ThreeVector GetGlobalExitNormal(const G4ThreeVector& point,
                                            G4bool* valid) override;

    /**
     * Computes the isotropic safety for 'Target'.
     *  @param[in] pGlobalpoint Point in global coordinates system.
     *  @returns The isotropic safety value.
     */
    G4double TargetSafetyFromPoint( const G4ThreeVector& pGlobalpoint );
};

#endif
