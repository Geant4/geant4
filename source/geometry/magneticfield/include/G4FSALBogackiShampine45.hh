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
// G4FSALBogackiShampine45
//
// Class description:
//
// Bogacki-Shampine - 8 - 5(4) FSAL stepper

// Author: Somnath Banerjee (CERN, Google Summer of Code 2015), 26.05.2015
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4FSAL_BOGACKI_SHAMPINE_45_HH
#define G4FSAL_BOGACKI_SHAMPINE_45_HH

#include "G4VFSALIntegrationStepper.hh"

/**
 * @brief G4FSALBogackiShampine45 is an integrator of particle's equation of
 * motion based on the Bogacki-Shampine - 8 - 5(4) FSAL implementation.
 */

class G4FSALBogackiShampine45 : public G4VFSALIntegrationStepper
{
  public:

    /**
     * Constructor for G4FSALBogackiShampine45.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     *  @param[in] numberOfVariables The number of integration variables.
     *  @param[in] primary Flag for initialisation of the auxiliary stepper.
     */
    G4FSALBogackiShampine45(G4EquationOfMotion* EqRhs,
                            G4int numberOfVariables = 6,
                            G4bool primary = true);

    /**
     * Destructor.
     */
    ~G4FSALBogackiShampine45() override;
    
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4FSALBogackiShampine45(const G4FSALBogackiShampine45&) = delete;
    G4FSALBogackiShampine45& operator=(const G4FSALBogackiShampine45&) = delete;

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
     *  @param[out] nextDydx Last derivatives array for the next step.
     */
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[],
                        G4double nextDydx[]) override ;
    
    /**
     * Calculates the output at the tau fraction of step.
     *  @param[in] yInput Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[out] yOut Interpolation output.
     *  @param[in] Step The given step size.
     *  @param[in] tau The tau fraction of the step.
     */
    void interpolate( const G4double yInput[],
                      const G4double dydx[],
                            G4double yOut[],
                            G4double Step,
                            G4double tau ) ;

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override;

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }
    
  private:
    
    /**
     * Init method used in constructor.
     */
    void PrepareConstants(); 
   
    /** Working arrays -- used during stepping. */
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9, *ak10, *ak11,
             *DyDx, *yTemp, *yIn;
    G4double *pseudoDydx_for_DistChord;

    G4double fLastStepLength = -1.0;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError; // for DistChord calculations

    /** Working array for interpolation. */
    G4double b[12];
   
    G4FSALBogackiShampine45* fAuxStepper = nullptr;

    static G4bool fPreparedConstants;
    static G4double bi[12][7];
};

#endif
