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
// G4CashKarpRKF45
//
// Class description:
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.

// Authors: J.Apostolakis, V.Grichine (CERN), 30.01.1997
// -------------------------------------------------------------------
#ifndef G4CASHKARP_RKF45_HH
#define G4CASHKARP_RKF45_HH

#include "G4MagIntegratorStepper.hh"

/**
 * @brief G4CashKarpRKF45 implements the Cash-Karp Runge-Kutta-Fehlberg
 * 4/5 method, an embedded fourth order method (giving fifth-order accuracy)
 * for the solution of an ODE. Two different fourth order estimates are
 * calculated; their difference gives an error estimate. 
 * It is used to integrate the equations of the motion of a particle 
 * in a magnetic field.
 */

class G4CashKarpRKF45 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4CashKarpRKF45.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4CashKarpRKF45( G4EquationOfMotion* EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary = true );

    /**
     * Destructor.
     */
    ~G4CashKarpRKF45() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4CashKarpRKF45(const G4CashKarpRKF45&) = delete;
    G4CashKarpRKF45& operator=(const G4CashKarpRKF45&) = delete;
 
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
     * Returns the distance from chord line.
     */
    G4double DistChord() const override; 

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kCashKarpRKF45".
     */
    inline G4StepperType StepperType() const override { return kCashKarpRKF45; }

  private:

    /** Scratch space. */
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;

    G4double fLastStepLength = 0.0;

    /** For DistChord calculations. */
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;

    G4CashKarpRKF45* fAuxStepper = nullptr;
};

#endif
