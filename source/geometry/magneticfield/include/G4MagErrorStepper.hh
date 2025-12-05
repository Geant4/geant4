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
// G4MagErrorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field.

// Author: W.Wander (MIT), 09.12.1997
// --------------------------------------------------------------------
#ifndef G4MAGERRORSTEPPER_HH
#define G4MAGERRORSTEPPER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

/**
 * @brief G4MagErrorStepper is an abstract base class for integrator of
 * particle's equation of motion, used in tracking in space dependent
 * magnetic field.
 */

class G4MagErrorStepper : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4MagErrorStepper.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] numberOfVariables The number of state variables.
     */
    G4MagErrorStepper(G4EquationOfMotion*EqRhs,
                      G4int numberOfVariables,
                      G4int numStateVariables = 12);

    /**
     * Destructor.
     */
    ~G4MagErrorStepper() override;
  
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4MagErrorStepper(const G4MagErrorStepper&) = delete;
    G4MagErrorStepper& operator=(const G4MagErrorStepper&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values y[0 to 6].
     * Outputs yout[] and its estimated error yerr[].
     *  @param[in] y Starting values array of integration variables.
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
     * without error calculation. To be implemented in concrete derived classes.
     */
    virtual void DumbStepper( const G4double y[],
                              const G4double dydx[],
                                    G4double h,
                                    G4double yout[] ) = 0;

    /**
     * Estimates the maximum distance of curved solution and chord.
     */
    G4double DistChord() const override;

  private:

    /** Data stored in order to find the chord. */
    G4ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
 
    /** Arrays used only for temporary storage; they are allocated at the
        class level only for efficiency, so that calls to new and delete are
        not made in Stepper(). */
    G4double *yInitial, *yMiddle, *dydxMid, *yOneStep;
};

#include  "G4MagErrorStepper.icc"

#endif
