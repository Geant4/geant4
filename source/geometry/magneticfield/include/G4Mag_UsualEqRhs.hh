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
// G4Mag_UsualEqRhs
//
// Class description:
//
// This is the standard right-hand side for equation of motion.
// The only case another is required is when using a moving reference
// frame ... or extending the class to include additional Forces,
// eg an electric field

// Author: John Apostolakis (CERN), 13.01.1997
// --------------------------------------------------------------------
#ifndef G4MAG_USUAL_EQRHS
#define G4MAG_USUAL_EQRHS

#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"

class G4MagneticField;

/**
 * @brief G4Mag_UsualEqRhs defines the standard right-hand side
 * for equation of motion.
 */

class G4Mag_UsualEqRhs : public G4Mag_EqRhs
{
  public:

    /**
     * Constructor for G4Mag_UsualEqRhs.
     *  @param[in] MagField Pointer to the associated magnetic field.
     */
    G4Mag_UsualEqRhs( G4MagneticField* MagField );

    /**
     * Default Destructor.
     */
    ~G4Mag_UsualEqRhs() override = default;
      // Constructor and destructor. No actions.

    /**
     * Calculates the value of the derivative, given the value of the field.
     *  @param[in] y Coefficients array.
     *  @param[in] B Field value.
     *  @param[out] dydx Derivatives array.
     */
    void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[3],
                                  G4double dydx[] ) const override;

    /**
     * Sets the charge momentum mass value.
     */
    void SetChargeMomentumMass( G4ChargeState particleCharge,
                                G4double MomentumXc,
                                G4double mass ) override;

    /**
     * Returns the equation of motion type ID, i.e. "kEqMagnetic".
     */
    inline G4EquationType GetEquationType() const override { return kEqMagnetic; }
};

#endif
