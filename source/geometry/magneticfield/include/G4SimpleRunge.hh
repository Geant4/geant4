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
// G4SimpleRunge
//
// Class description:
//
// Simple Runge:
//
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
// Second order solver.
// Takes the derivative at a position to be assumed at the middle of the
// Step and adds it to the current position.

// Author: W.Wander (MIT), 12.09.1997
// -------------------------------------------------------------------
#ifndef G4SIMPLERUNGE_HH
#define G4SIMPLERUNGE_HH

#include "G4MagErrorStepper.hh"

/**
 * @brief G4SimpleRunge implements a simple Runge stepper for magnetic field
 * with 2nd order solver.
 */

class G4SimpleRunge : public G4MagErrorStepper
{
  public:

    /**
     * Constructor for G4SimpleRunge.
     *  @param[in] EquationRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     */
    G4SimpleRunge(G4EquationOfMotion* EquationRhs,
                  G4int numberOfVariables = 6);

    /**
     * Destructor.
     */
    ~G4SimpleRunge() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4SimpleRunge(const G4SimpleRunge&) = delete;
    G4SimpleRunge& operator=(const G4SimpleRunge&) = delete;

    /**
     * The stepper for the Runge Kutta integration, but performing a 'dump' step
     * without error calculation.
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
     * Returns the order, 2, of integration.
     */
    inline G4int IntegratorOrder() const override { return 2; }

    /**
     * Returns the stepper type-ID, "kSimpleRunge".
     */
    inline G4StepperType StepperType() const override { return kSimpleRunge; }

  private:

    G4int fNumberOfVariables = 0;

    G4double* dydxTemp = nullptr;
    G4double* yTemp = nullptr;
};

#endif
