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
// G4BorisDriver
//
// Class description:
//
// G4BorisDriver is a driver class using the second order Boris 
// method to integrate the equation of motion.

// Author: Divyansh Tiwari (CERN, Google Summer of Code 2022), 05.11.2022
// Supervision: John Apostolakis (CERN), Renee Fatemi, Soon Yung Jun (FNAL)
// --------------------------------------------------------------------
#ifndef G4BORIS_DRIVER_HH
#define G4BORIS_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4BorisScheme.hh"
#include "G4ChordFinderDelegate.hh"

/**
 * @brief G4BorisDriver is a driver class using the second order Boris 
 * method to integrate the equation of motion.
 */

class G4BorisDriver : public G4VIntegrationDriver,
                      public G4ChordFinderDelegate<G4BorisDriver>
{
  public:

    /**
     * Constructor for G4BorisDriver.
     *  @param[in] hminimum The minumum allowed step.
     *  @param[in] Boris Pointer to the Boris motion algorithm.
     *  @param[in] numberOfComponents The number of integration variables.
     *  @param[in] verbosity Flag for verbosity.
     */
    G4BorisDriver( G4double hminimum,
                   G4BorisScheme* Boris,
                   G4int numberOfComponents = 6,
                   G4bool verbosity = false);
   
    /**
     * Default Destructor.
     */
    ~G4BorisDriver() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4BorisDriver(const G4BorisDriver&) = delete;
    G4BorisDriver& operator=(const G4BorisDriver&) = delete;

    // 1. Core methods that advance the integration

    /**
     * Advances integration accurately by relative accuracy better than 'eps'.
     *  @param[in,out] track The current track in field.
     *  @param[in] stepLen Proposed step length.
     *  @param[in] epsilon Requested accuracy, y_err/hstep.
     *  @param[in] beginStep Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    G4bool AccurateAdvance(G4FieldTrack& track,
                           G4double stepLen,
                           G4double eps,
                           G4double beginStep = 0) override;

    /**
     * Attempts one integration step, and returns estimated error 'dyerr'.
     * It does not ensure accuracy.
     *  @param[in,out] y_val The current track in field.
     *  @param[in] dydx dydx array.
     *  @param[in] hstep Proposed step length.
     *  @param[out] missDist Estimated sagitta distance.
     *  @param[out] dyerr Estimated error.
     *  @returns true if integration succeeds.
     */
    G4bool QuickAdvance(G4FieldTrack& y_val,   // In/Out
                        const G4double dydx[],    
                        G4double hstep,
                        G4double& missDist,    // Out: estimated sagitta
                        G4double& dyerr) override;
    /**
     * Takes one Step that is as large as possible while satisfying the
     * accuracy criterion.
     *  @param[in,out] yCurrentState The current track state, y.
     *  @param[in,out] curveLength Step start, x.
     *  @param[in] htry Step to attempt.
     *  @param[in] epsilon_rel The relative accuracy.
     *  @param[in] restMass Mass value for computing velocity.
     *  @param[in] charge Charge value for computing momentum.
     *  @param[out] hdid Step achieved.
     *  @param[out] hnext Proposed next step.
     *  @returns true if integration succeeds.
     */
    void OneGoodStep(G4double yCurrentState[],  // In/Out: state ('y')
                     G4double& curveLength,     // In/Out: 'x'
                     G4double htry,             // step to attempt
                     G4double epsilon_rel,      // relative accuracy
                     G4double restMass,         
                     G4double charge,
                     G4double& hdid,            // Out: step achieved
                     G4double& hnext);          // Out: proposed next step
   
    // 2. Methods needed to co-work with G4ChordFinder

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    inline G4double AdvanceChordLimited(G4FieldTrack& track,
                                        G4double hstep,
                                        G4double eps,
                                        G4double chordDistance) override;
    /**
     * Dispatch interface method for initialisation/reset of driver.
     */
    inline void OnStartTracking() override;

    /**
     * Dispatch interface method for computing step. Does nothing here.
     */
    inline void OnComputeStep(const G4FieldTrack*) override;

    // 3. Does the method redo integrations when called to obtain values for
    //    internal, smaller intervals? (when needed to identify an intersection)

    /**
     * The driver implements re-integration. Returns true.
     * It would be false if it just used interpolation to provide a result.
     */
    inline G4bool DoesReIntegrate() const override;

    // 4. Relevant for calculating a new step size to achieve required accuracy

    /**
     * Computes a step size for the next step, taking the last step's
     * normalised error 'errMaxNorm'.
     *  @param[in] errMaxNorm The normalised error on last step.
     *  @param[in] hstepCurrent The current proposed step.
     *  @returns The step size for the next step.
     */
    inline G4double ComputeNewStepSize(G4double errMaxNorm,
                                       G4double hstepCurrent) override;

    /**
     * Methods to calculate the next step size given the square of the
     * relative error.
     */
    G4double ShrinkStepSize2(G4double h, G4double error2) const;
    G4double GrowStepSize2(G4double h, G4double error2) const;
   
    // 5. Auxiliary Methods ...

    /**
     * Getters for derivatives.
     */
    void GetDerivatives( const G4FieldTrack& track,
                               G4double dydx[] ) const override;
    void GetDerivatives( const G4FieldTrack& track,
                               G4double dydx[],
                               G4double field[] ) const override;

    /**
     * Setter and getter for verbosity.
     */
    inline void SetVerboseLevel(G4int level) override;
    inline G4int GetVerboseLevel() const override;

    /**
     * Getters for the equation of motion.
     */
    inline G4EquationOfMotion* GetEquationOfMotion() override;
    inline const G4EquationOfMotion* GetEquationOfMotion() const;

    /**
     * Setter for the equation of motion. Issues an exception, as not
     * foreseen to change equation of motion for the Boris stepper.
     */
    void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    void StreamInfo( std::ostream& os ) const override;
   
    /**
     * Accessors for stepper. Not relevant for Boris and other non-RK methods.
     */
    inline const G4MagIntegratorStepper* GetStepper() const override;
    inline G4MagIntegratorStepper* GetStepper() override;

  private:

    inline G4int GetNumberOfVariables() const;

    inline void CheckStep(const G4ThreeVector& posIn,                              
                          const G4ThreeVector& posOut,
                                G4double hdid) const;
   
  private:

    // INVARIANTS -- remain unchanged during tracking / integration 
    // Parameters
    G4double fMinimumStep;
    G4bool   fVerbosity;

    // State -- The core stepping algorithm
    G4BorisScheme* boris;
    
    // STATE -- intermediate state (to avoid creation / churn )
    G4double yIn[G4FieldTrack::ncompSVEC],
             yMid[G4FieldTrack::ncompSVEC],
             yOut[G4FieldTrack::ncompSVEC],
             yError[G4FieldTrack::ncompSVEC];

    G4double yCurrent[G4FieldTrack::ncompSVEC];

    // - Unused 2022.11.03:   
    // G4double derivs[2][6][G4FieldTrack::ncompSVEC];
    // const G4int interval_sequence[2];
   
    // INVARIANTS -- Parameters for ensuring that one call has finite number of integration steps
    static constexpr G4int    fMaxNoSteps = 300; 
    static constexpr G4double fSmallestFraction= 1e-12; // To avoid FP underflow !  ( 1.e-6 for single prec)

    static constexpr G4int    fIntegratorOrder= 2; //  2nd order method -- needed for error control
    static constexpr G4double fSafetyFactor = 0.9; //

    static constexpr G4double fMaxSteppingIncrease= 10.0; //  Increase no more than 10x   
    static constexpr G4double fMaxSteppingDecrease= 0.1;  //  Reduce   no more than 10x
    static constexpr G4double fPowerShrink = -1.0 / fIntegratorOrder;
    static constexpr G4double fPowerGrow   = -1.0 / (1.0 + fIntegratorOrder);

    static const G4double fErrorConstraintShrink;
    static const G4double fErrorConstraintGrow;
   
    using ChordFinderDelegate = G4ChordFinderDelegate<G4BorisDriver>;
};

#include "G4BorisDriver.icc"

#endif
