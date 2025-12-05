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
// G4DormandPrinceRK78
//
// Class description:
//
// Dormand-Prince 8(7)13M non-FSAL RK method, a 13 stage embedded
// explicit Runge-Kutta method, using a pair of 7th and 8th order formulae.
//
// Paper proposing this RK scheme:
// P.J. Prince, J.R. Dormand, "High order embedded Runge-Kutta formulae",
// Journal of Computational and Applied Mathematics, Volume 7, Issue 1, 1981,
// Pages 67-75, ISSN 0377-0427, DOI: 10.1016/0771-050X(81)90010-3

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 28.06.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4DORMAND_PRINCE_RK78_HH
#define G4DORMAND_PRINCE_RK78_HH

#include "G4MagIntegratorStepper.hh"

/**
 * @brief G4DormandPrinceRK78 implements the Dormand-Prince 8(7)13M non-FSAL
 * Runge-Kutta method, a 13 stage embedded explicit Runge-Kutta method, using
 * a pair of 7th and 8th order formulae.
 */

class G4DormandPrinceRK78 : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4DormandPrince745.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4DormandPrinceRK78(G4EquationOfMotion* EqRhs,
                        G4int numberOfVariables = 6,
                        G4bool primary = true);

    /**
     * Destructor.
     */
    ~G4DormandPrinceRK78() override;
    
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4DormandPrinceRK78(const G4DormandPrinceRK78&) = delete;
    G4DormandPrinceRK78& operator=(const G4DormandPrinceRK78&) = delete;

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
                        G4double yerr[]) override ;

    /**
     * Returns the distance from chord line.
     */
    G4double  DistChord() const override;

    /**
     * Returns the order, 7, of integration.
     */
    inline G4int IntegratorOrder() const override { return 7; }

    /**
     * Returns the stepper type-ID, "kDormandPrinceRK78".
     */
    inline G4StepperType StepperType() const override { return kDormandPrinceRK78; }
    
  private :  
    
    G4double *ak2,   *ak3,  *ak4,  *ak5,  *ak6,  *ak7, *ak8,
             *ak9,   *ak10, *ak11, *ak12, *ak13,
             *yTemp, *yIn;

    G4double fLastStepLength = -1.0;

    /** For DistChord() calculations. */
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;

    G4DormandPrinceRK78* fAuxStepper = nullptr;
};

#endif
