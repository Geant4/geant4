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
// G4ErrorMag_UsualEqRhs
//
// Class description:
//
// Serves to reverse the magnetic field when propagation is backwards
// for error propagation.

// Author: Pedro Arce (CIEMAT), September 2004.
// --------------------------------------------------------------------
#ifndef G4ERRORMAG_USUALEQRHS_HH
#define G4ERRORMAG_USUALEQRHS_HH

#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"

/**
 * @brief G4ErrorMag_UsualEqRhs serves to reverse the magnetic field when
 * propagation is backwards. It is used for error propagation.
 */

class G4ErrorMag_UsualEqRhs : public G4Mag_UsualEqRhs
{
  public:

    /**
     * Constructor for G4ErrorMag_UsualEqRhs.
     *  @param[in] MagField Pointer to the magnetic field.
     */
    G4ErrorMag_UsualEqRhs( G4MagneticField* MagField );

    /**
     * Default Destructor.
     */
    ~G4ErrorMag_UsualEqRhs() override = default;

    /**
     * Calculates the value of the derivative, given the value of the
     * magnetic field. Reverses dedx if propagation is backwards.
     *  @param[in] y Coefficients array.
     *  @param[in] B Field value.
     *  @param[out] dydx Derivatives array.
     */
    void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[3],
                                  G4double dydx[] ) const override;
};

#endif
