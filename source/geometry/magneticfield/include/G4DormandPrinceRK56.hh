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
// G4DormandPrinceRK56
//
// Class description:
//
// Dormand-Prince RK 6(5) non-FSAL method

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 26.06.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4DORMAND_PRINCE_RK56_HH
#define G4DORMAND_PRINCE_RK56_HH

#include "G4MagIntegratorStepper.hh"

/**
 * @brief G4DormandPrinceRK56 implements the 6(5) embedded Runge-Kutta
 * non-FSAL method.
 */

class G4DormandPrinceRK56 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4DormandPrinceRK56.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4DormandPrinceRK56( G4EquationOfMotion* EqRhs,
                         G4int numberOfVariables = 6,
                         G4bool primary = true ) ;
    
    /**
     * Destructor.
     */
    ~G4DormandPrinceRK56() override ;
    
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4DormandPrinceRK56(const G4DormandPrinceRK56&) = delete;
    G4DormandPrinceRK56& operator=(const G4DormandPrinceRK56&) = delete;
    
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
                        G4double yerr[] ) override ;
    
    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 5, of integration.
     */
    inline G4int IntegratorOrder() const override { return 5; }

    /**
     * Returns the stepper type-ID, "kDormandPrinceRK56".
     */
    inline G4StepperType StepperType() const override { return kDormandPrinceRK56; }
    
    /**
     * Prepares the interpolant and calculates the extra stages.
     * Fifth order interpolant with one extra function evaluation per step.
     */
    void SetupInterpolate_low( const G4double yInput[],
                               const G4double dydx[],
                               const G4double Step );

    /**
     * Wrappers for SetupInterpolate_low() above.
     */
    inline void SetupInterpolate( const G4double yInput[],
                                  const G4double dydx[],
                                  const G4double Step )
    {
      SetupInterpolate_low( yInput, dydx, Step);
    }
    inline void SetupInterpolation()
    {
      SetupInterpolate( fLastInitialVector, fLastDyDx, fLastStepLength);
    }

    /**
     * Calculates the output at the tau fraction of Step.
     */
    void Interpolate_low( const G4double yInput[],
                          const G4double dydx[],
                          const G4double Step,
                                G4double yOut[],
                                G4double tau );

    /**
     * Wrappers for Interpolate_low() above.
     */
    inline void Interpolate( const G4double yInput[],
                             const G4double dydx[],
                             const G4double Step,
                                   G4double yOut[],
                                   G4double tau )
    {
      Interpolate_low( yInput, dydx, Step, yOut, tau);
    }
    inline void Interpolate( G4double tau, G4double yOut[])
    {
      Interpolate( fLastInitialVector, fLastDyDx, fLastStepLength, yOut, tau );
    }

    /**
     * Prepares the interpolant and calculates the extra stages.
     * Sixth order interpolant with 3 additional stages per step.
     */
    void SetupInterpolate_high( const G4double yInput[],
                                const G4double dydx[],
                                const G4double Step );
    
    /**
     * Calculates the output at the tau fraction of Step, using
     * the polynomial coefficients and the respective stages.
     */
    void Interpolate_high( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step,
                                 G4double yOut[],
                                 G4double tau );
    
  private:

    /** For storing intermediate 'k' values in stepper. */
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9;

    /** For the additional stages of Interpolant. */
    G4double *ak10_low, *ak10, *ak11, * ak12;

    G4double *yTemp, *yIn;
    
    G4double fLastStepLength = -1.0;

    /** For DistChord() calculations. */
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
    
    G4DormandPrinceRK56* fAuxStepper = nullptr;
};

#endif
