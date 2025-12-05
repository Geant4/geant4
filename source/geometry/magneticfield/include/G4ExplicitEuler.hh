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
// G4ExplicitEuler
//
// Class description:
//
// Explicit Euler stepper for magnetic field: x_1 = x_0 + h * dx_0.
// The most simple approach for solving linear differential equations.
// Take the current derivative and add it to the current position.

// Author: W.Wander (MIT), 12.09.1997
// -------------------------------------------------------------------
#ifndef G4EXPLICITEULER_HH
#define G4EXPLICITEULER_HH

#include "G4MagErrorStepper.hh"

/**
 * @brief G4ExplicitEuler implements an Explicit Euler stepper for magnetic
 * field: x_1 = x_0 + h * dx_0. The most simple approach for solving linear
 * differential equations. Takes the current derivative and adds it to the
 * current position.
 */

class G4ExplicitEuler : public G4MagErrorStepper
{
  public:

    /**
     * Constructor for G4ExplicitEuler.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     */
    G4ExplicitEuler(G4EquationOfMotion* EqRhs,
                    G4int numberOfVariables = 6) ;

    /**
     * Default Destructor.
     */
    ~G4ExplicitEuler() override = default;

    /**
     * The stepper function for the integration.
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     */
    void DumbStepper( const G4double y[],
                       const G4double dydx[],
                             G4double h,
                             G4double yout[] ) override;

    /**
     * Returns the order, 1, of integration.
     */
    inline G4int IntegratorOrder() const override { return 1; }

    /**
     * Returns the stepper type-ID, "kExplicitEuler".
     */
    inline G4StepperType StepperType() const override { return kExplicitEuler; }
};

#endif
