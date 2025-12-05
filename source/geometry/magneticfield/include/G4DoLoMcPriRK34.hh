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
// G4DoLoMcPriRK34
//
// Class description:
//
// Dormand-Lockyer-McGorrigan-Prince-6-3-4 non-FSAL method
// ( 6 stage, 3rd & 4th order embedded RK method )

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 07.07.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef DOLO_MCPRI_RK34_HH
#define DOLO_MCPRI_RK34_HH

#include "G4MagIntegratorStepper.hh"

/**
 * @brief G4DoLoMcPriRK34 implements the Dormand-Lockyer-McGorrigan-Prince-6-3-4
 * non-FSAL method ( 6 stage, 3rd & 4th order embedded Runge-Kutta method ).
 */

class G4DoLoMcPriRK34 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4DoLoMcPriRK34.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4DoLoMcPriRK34( G4EquationOfMotion* EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary = true );

    /**
     * Destructor.
     */
    ~G4DoLoMcPriRK34() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4DoLoMcPriRK34(const G4DoLoMcPriRK34&) = delete;
    G4DoLoMcPriRK34& operator=(const G4DoLoMcPriRK34&) = delete; 

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
     * Interface method for interpolation setup. Does nothing here.
     */
    inline void SetupInterpolation() {}

    /**
     * Calculates the output at the tau fraction of Step.
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] Step The given step size.
     *  @param[out] yOut Interpolation output.
     *  @param[out] tau Fraction of step.
     */
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                            G4double yOut[],
                            G4double tau );
    void Interpolate( G4double tau,
                      G4double yOut[]);
    
    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 3, of integration.
     */
    inline G4int IntegratorOrder() const override { return 3; }


    /**
     * Returns the stepper type-ID, "kDoLoMcPriRK34".
     */
    inline G4StepperType StepperType() const override { return kDoLoMcPriRK34; }
    
  private :
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;
    
    G4double fLastStepLength = -1.0;

    /** For DistChord calculations. */
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
    
    G4DoLoMcPriRK34* fAuxStepper = nullptr;
};

#endif
