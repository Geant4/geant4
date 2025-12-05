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
// G4MagInt_Driver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.

// Author: Vladimir Grichine (CERN), 07.10.1996 - Created
//         W.Wander (MIT), 28.01.1998 - Added ability for low order integrators
// --------------------------------------------------------------------
#ifndef G4MAGINT_DRIVER_HH
#define G4MAGINT_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinderDelegate.hh"

/**
 * @brief G4MagInt_Driver provides a driver that talks to the Integrator
 * Stepper and insures that the error is within acceptable bounds.
 */

class G4MagInt_Driver : public G4VIntegrationDriver,
                        public G4ChordFinderDelegate<G4MagInt_Driver>
{
  public:

    /**
     * Constructor for G4MagInt_Driver.
     *  @param[in] hminimum The minumum allowed step.
     *  @param[in] pItsStepper Pointer to the integrator stepper.
     *  @param[in] numberOfComponents The number of integration variables.
     *  @param[in] statisticsVerbosity Flag for verbosity.
     */
    G4MagInt_Driver(G4double hminimum,
                    G4MagIntegratorStepper* pItsStepper,
                    G4int numberOfComponents = 6,
                    G4int statisticsVerbosity = 0);

    /**
     * Destructor. Provides statistics if verbosity level is greater than 1.
     */
    ~G4MagInt_Driver() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4MagInt_Driver(const G4MagInt_Driver&) = delete;
    G4MagInt_Driver& operator=(const G4MagInt_Driver&) = delete;

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] stepMax Proposed maximum step length.
     *  @param[in] epsStep Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    inline G4double AdvanceChordLimited(G4FieldTrack& track, 
                                        G4double stepMax, 
                                        G4double epsStep,
                                        G4double chordDistance) override;

    /**
     * Dispatch interface method for initialisation/reset of driver.
     */
    inline void OnStartTracking() override;

    /**
     * Dispatch interface method for computing step. Does nothing here.
     */
    inline void OnComputeStep(const G4FieldTrack* = nullptr) override {}

    /**
     * The driver implements re-integration, so returns true.
     */
    G4bool DoesReIntegrate() const override { return true; }
   
    /**
     * Advances integration accurately by relative accuracy better than 'eps'.
     *  @param[in,out] y_current The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] hinitial Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    G4bool AccurateAdvance(G4FieldTrack& y_current,
                           G4double hstep,
                           G4double eps,  // Requested y_err/hstep
                           G4double hinitial = 0.0) override;

    /**
     * Attempts one integration step, and returns estimated error 'dyerr'.
     * It does not ensure accuracy.
     *  @param[in,out] y_val The current track in field.
     *  @param[in] dydx dydx array.
     *  @param[in] hstep Proposed step length.
     *  @param[out] dchord_step Estimated sagitta distance.
     *  @param[out] dyerr Estimated error.
     *  @returns true if integration succeeds.
     */
    G4bool QuickAdvance(G4FieldTrack& y_val,   // In/Out
                        const G4double dydx[],
                        G4double hstep,
                        G4double& dchord_step,
                        G4double& dyerr) override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    void StreamInfo( std::ostream& os ) const override;

    /**
     * Attempts one integration step, and returns estimated error 'dyerr'.
     * It does not ensure accuracy.
     *  @param[in,out] y_posvel The current track in field.
     *  @param[in] dydx dydx array.
     *  @param[in] hstep Proposed step length.
     *  @param[out] dchord_step Estimated sagitta distance.
     *  @param[out] dyerr_pos_sq Estimated error in position.
     *  @param[out] dyerr_mom_rel_sq Estimated error in momentum
     *              (normalised: Delta_Integration(p^2)/(p^2)).
     *  @returns true if integration succeeds.
     */
    G4bool QuickAdvance(G4FieldTrack& y_posvel,   // In/Out
                        const G4double dydx[],
                        G4double hstep,           // In
                        G4double& dchord_step,
                        G4double& dyerr_pos_sq,
                        G4double& dyerr_mom_rel_sq );

    /**
     * Accessors.
     */
    inline G4double GetHmin() const;
    inline G4double Hmin() const;     // Obsolete
    inline G4double GetSafety() const;
    inline G4double GetPshrnk() const;
    inline G4double GetPgrow() const;
    inline G4double GetErrcon() const;
    void GetDerivatives(const G4FieldTrack& y_curr,            // INput
                              G4double dydx[]) const override; // OUTput
    void GetDerivatives(const G4FieldTrack& track,
                              G4double dydx[],
                              G4double field[]) const override;

    /**
     * Getter and setter for the equation of motion.
     */
    G4EquationOfMotion* GetEquationOfMotion() override;
    void SetEquationOfMotion(G4EquationOfMotion* equation) override;
   
    /**
     * Sets a new stepper 'pItsStepper' for this driver. Then it calls
     * ResetParameters() to update its parameters accordingly.
     */
    void RenewStepperAndAdjust(G4MagIntegratorStepper* pItsStepper) override;

    /**
     * Resets the qarameters according to the new provided safety value.
     *  i) sets the exponents (pgrow & pshrnk), using the current order;
     * ii) sets the safety and calculates "errcon" according to the above values.
     */
    inline void ReSetParameters(G4double new_safety = 0.9);

    /**
     * Modifiers. When setting safety or pgrow, errcon will be set
     * to a compatible value.
     */
    inline void SetSafety(G4double valS);
    inline void SetPshrnk(G4double valPs);
    inline void SetPgrow (G4double valPg);
    inline void SetErrcon(G4double valEc);
    inline G4double ComputeAndSetErrcon();

    /**
     * Accessors for the integrator stepper.
     */
    const G4MagIntegratorStepper* GetStepper() const override;
    G4MagIntegratorStepper* GetStepper() override;

    /**
     * Takes one Step that is as large as possible while satisfying the
     * accuracy criterion of: yerr < eps * |y_end-y_start|.
     *  @param[in,out] ystart The current track state, y.
     *  @param[in] dydx The derivatives array.
     *  @param[in,out] x Step start, x.
     *  @param[in] htry Step to attempt.
     *  @param[in] eps The relative accuracy.
     *  @param[out] hdid Step achieved.
     *  @param[out] hnext Proposed next step.
     *  @returns true if integration succeeds.
     */
    void OneGoodStep(G4double ystart[], // Like old RKF45step()
                     const G4double dydx[],
                     G4double& x,
                     G4double htry,
                     G4double eps,
                     G4double& hdid,
                     G4double& hnext ) ;

    /**
     * Takes the last step's normalised error and calculates a step size
     * for the next step. Does it limit the next step's size within a factor
     * of the current?
     * --  DOES NOT limit for very bad steps
     * --  DOES     limit for very good (x5).
     */
    G4double ComputeNewStepSize(G4double errMaxNorm, // normalised
                                G4double hstepCurrent) override;

    /**
     * Taking the last step's normalised error, calculates a step size for
     * the next step. Does not limit the next step's size within a factor of
     * the current one when *reducing* the size, i.e. for badly failing steps.
     */
    G4double ComputeNewStepSize_WithoutReductionLimit(G4double errMaxNorm,
                                                      G4double hstepCurrent);
                                                      
    /**
     * Taking the last step's normalised error, calculates a step size for
     * the next step. Limits the next step's size within a range around the
     * current one.
     */
    G4double ComputeNewStepSize_WithinLimits(G4double errMaxNorm, // normalised
                                             G4double hstepCurrent);

    /**
     * Modifier and accessor for the maximum number of steps that can be taken
     * for the integration of a single segment, i.e. a single call to
     * AccurateAdvance().
     */
    inline G4int GetMaxNoSteps() const;
    inline void SetMaxNoSteps(G4int val);

    /**
     * More modifiers and accessors.
     */
    inline void SetHmin(G4double newval);
    void SetVerboseLevel(G4int newLevel) override;
    G4int GetVerboseLevel() const override;
    inline G4double GetSmallestFraction() const;
    void SetSmallestFraction( G4double val );

  protected:

    /**
     * Loggers, issuing warnings for undesirable situations.
     */
    void WarnSmallStepSize(G4double hnext, G4double hstep,
                           G4double h, G4double xDone,
                           G4int noSteps);
    void WarnTooManyStep(G4double x1start, G4double x2end, G4double xCurrent);
    void WarnEndPointTooFar(G4double endPointDist,
                            G4double hStepSize ,
                            G4double epsilonRelative,
                            G4int debugFlag);

    /**
     * Loggers for verbosity printouts.
     */
    void PrintStatus(const G4double* StartArr,
                           G4double xstart,
                     const G4double* CurrentArr,
                           G4double xcurrent,
                           G4double requestStep,
                           G4int subStepNo);
    void PrintStatus(const G4FieldTrack& StartFT,
                     const G4FieldTrack& CurrentFT,
                           G4double requestStep,
                           G4int subStepNo);
    void PrintStat_Aux(const G4FieldTrack& aFieldTrack,
                             G4double requestStep,
                             G4double actualStep,
                             G4int subStepNo,
                             G4double subStepSize,
                             G4double dotVelocities);
    /**
     * Reports on the number of steps, maximum errors etc.
     */
    void PrintStatisticsReport();

#ifdef QUICK_ADV_TWO
    G4bool QuickAdvance(      G4double  yarrin[],     // In
                        const G4double  dydx[],  
                              G4double  hstep,        
                              G4double  yarrout[],    // Out
                              G4double& dchord_step,  // Out
                              G4double& dyerr );      // in length
#endif

  private:

    // ---------------------------------------------------------------
    //  INVARIANTS 

    /** Minimum Step allowed in a Step (in absolute units). */
    G4double fMinimumStep = 0.0;

    /** Smallest fraction of (existing) curve length, in relative units.
        Below this fraction the current step will be the last. */
    G4double fSmallestFraction = 1.0e-12;    // Expected range 1e-12 to 5e-15

    /** Variables in integration. */
    const G4int fNoIntegrationVariables = 0;

    /** Minimum number for FieldTrack. */
    const G4int fMinNoVars = 12;

    /** Full number of variable. */
    const G4int fNoVars = 0;

    /** Default maximum number of steps is Base divided by the order of Stepper. */
    G4int fMaxNoSteps;
    G4int fMaxStepBase = 250;  // was 5000

    /** Parameters used to grow and shrink trial stepsize. */
    G4double safety;
    G4double pshrnk;   //  exponent for shrinking
    G4double pgrow;    //  exponent for growth
    G4double errcon;

    G4int fStatisticsVerboseLevel = 0;

    // ---------------------------------------------------------------
    // DEPENDENT Objects

    G4MagIntegratorStepper* pIntStepper = nullptr;

    // ---------------------------------------------------------------
    //  STATE

    /** Step Statistics. */
    unsigned long fNoTotalSteps=0, fNoBadSteps=0;
    unsigned long fNoSmallSteps=0, fNoInitialSmallSteps=0, fNoCalls=0;
    G4double fDyerr_max=0.0, fDyerr_mx2=0.0;
    G4double fDyerrPos_smTot=0.0, fDyerrPos_lgTot=0.0, fDyerrVel_lgTot=0.0;
    G4double fSumH_sm=0.0, fSumH_lg=0.0;

    /** Could be varied during tracking - to help identify issues. */
    G4int fVerboseLevel = 0;   // Verbosity level for printing (debug, ..)

    using ChordFinderDelegate = G4ChordFinderDelegate<G4MagInt_Driver>;
};

#include "G4MagIntegratorDriver.icc"

#endif
