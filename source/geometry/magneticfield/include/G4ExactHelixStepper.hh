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
// G4ExactHelixStepper
//
// Class description:
//
// Concrete class for particle motion in constant magnetic field.
// Helix a-la-Explicity Euler: x_1 = x_0 + helix(h)
// with helix(h) being a helix piece of length h.
// simplest approach for solving linear differential equations.
// Take the current derivative and add it to the current position.
//
// As the field is assumed constant, an error is not calculated.

// Author: John Apostolakis (CERN), 28.01.2005.
//         Implementation adapted from ExplicitEuler by W.Wander 
// --------------------------------------------------------------------
#ifndef G4EXACTHELIXSTEPPER_HH
#define G4EXACTHELIXSTEPPER_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4MagHelicalStepper.hh"
#include "G4Mag_EqRhs.hh"

/**
 * @brief G4ExactHelixStepper is a concrete class for particle motion in
 * constant magnetic field. Helix a-la-Explicity Euler: x_1 = x_0 + helix(h)
 * with helix(h) being a helix piece of length h.
 * As the field is assumed constant, an error is not calculated.
 */

class G4ExactHelixStepper : public G4MagHelicalStepper
{
  public:

    /**
     * Constructor for G4ExactHelixStepper.
     *  @param[in] EqRhs Pointer to the standard equation of motion.
     */
    G4ExactHelixStepper(G4Mag_EqRhs* EqRhs);

    /**
     * Default Destructor.
     */
    ~G4ExactHelixStepper() override = default;
  
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4ExactHelixStepper(const G4ExactHelixStepper&) = delete;
    G4ExactHelixStepper& operator=(const G4ExactHelixStepper&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Provides helix starting values y[0 to 6].
     * Outputs yout[] and ZERO estimated error yerr[]=0.
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yerr The estimated error.
     */
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) override;
  
    /**
     * Same as Stepper() function above, but should perform a 'dump' step
     * without error calculation. Assuming a constant field, the solution is
     * a helix. Should NOT be called; issues a fatal exception as the Stepper
     * must do all the work.
     */
    void DumbStepper( const G4double y[],
                            G4ThreeVector Bfld,
                            G4double h,
                            G4double yout[] ) override;
  
    /**
     * Estimates the maximum distance of curved solution and chord.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 1, of integration.
     */
    inline G4int IntegratorOrder() const override { return 1; }

    /**
     * Returns the stepper type-ID, "kExactHelixStepper".
     */
    inline G4StepperType StepperType() const override { return kExactHelixStepper; }

  private:

    /** Initial value of field at last step. */
    G4ThreeVector fBfieldValue;
};

#endif
