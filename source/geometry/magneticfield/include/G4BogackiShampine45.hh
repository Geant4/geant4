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
// G4BogackiShampine45
//
// Class description:
//
// An implementation of the embedded RK method from the following paper 
// by P. Bogacki and L. F. Shampine:
//    "An efficient Runge-Kutta (4,5) pair"
//    Comput. Math. with Appl., vol. 32, no. 6, pp. 15-28, Sep. 1996.
//
// An interpolation method provides the value of an intermediate
// point in a step -- if a step was sucessful.
//
// This version can provide the FSAL property of the method,
// which allows the reuse of the last derivative in the next step,
// but only by using the additional method GetLastDyDx() (an alternative
// interface for simpler use of FSAL is under development).

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 25.05.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef BOGACKI_SHAMPINE_45_HH
#define BOGACKI_SHAMPINE_45_HH

#include "G4MagIntegratorStepper.hh"

/**
 * @brief G4BogackiShampine45 is an integrator of particle's equation of
 * motion based on the Bogacki-Shampine method with FSAL property, allowing
 * the reuse of the last derivative in the next step.
 * This Stepper provides 'dense output'. After a successful step, it is
 * possible to obtain an estimate of the value of the function at an
 * intermediate point of the interval. This requires only two additional
 * evaluations of the derivative (and thus the field).
 */

class G4BogackiShampine45 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4BogackiShampine45.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4BogackiShampine45(G4EquationOfMotion* EqRhs,
                        G4int numberOfVariables = 6,
                        G4bool primary =  true);

    /**
     * Destructor.
     */
    ~G4BogackiShampine45() override;
    
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4BogackiShampine45(const G4BogackiShampine45&) = delete;
    G4BogackiShampine45& operator=(const G4BogackiShampine45&) = delete;

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
     * Setup all coefficients for interpolation.
     */
    void SetupInterpolationHigh();  // ( yInput, dydx, Step);
    inline void SetupInterpolation() { SetupInterpolationHigh(); }
    
    /**
     * Calculates the output at the tau fraction of step.
     *  @param[in] tau The tau fraction of the step.
     *  @param[out] yOut Interpolation output.
     */
    void InterpolateHigh( G4double tau, G4double yOut[] ) const;
    inline void Interpolate( G4double tau, 
                             G4double yOut[] ) { InterpolateHigh( tau, yOut); }

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kBogackiShampine45".
     */
    inline G4StepperType StepperType() const override { return kBogackiShampine45; }

    /**
     * Acccessor for dydx array.
     */
    void GetLastDydx( G4double dyDxLast[] );

    /**
     * Initialises the values of the bi[][] array.
     */
    void PrepareConstants();
   
  private:
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8,
             *ak9, *ak10, *ak11, *yTemp, *yIn;

    G4double *p[6];

    G4double fLastStepLength = -1.0;

    /** For DistChord() calculations. */
    G4double *fLastInitialVector, *fLastFinalVector, *fLastDyDx,
             *fMidVector, *fMidError;
    
    G4BogackiShampine45* fAuxStepper = nullptr;

    /** For chord - until interpolation is proven. */
    G4bool fPreparedInterpolation = false;

    /** Class constants. */
    static G4bool fPreparedConstants;   
    static G4double bi[12][7];
};

#endif
