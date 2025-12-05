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
// G4BulirschStoer driver
//
// Class description:
//
// G4IntegrationDriver<G4BulirschStoer> is a concrete driver class using
// Bulirsch-Stoer method to integrate the equation of motion.

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2016), 13.02.2018
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4BULIRSCH_STOER_DRIVER_HH
#define G4BULIRSCH_STOER_DRIVER_HH

#include "G4IntegrationDriver.hh"
#include "G4BulirschStoer.hh"
#include "G4ChordFinderDelegate.hh"

/**
 * @brief G4IntegrationDriver<G4BulirschStoer> is a concrete driver class
 * using the Bulirsch-Stoer method to integrate the equation of motion.
 */

template <>
class G4IntegrationDriver<G4BulirschStoer>: 
    public G4VIntegrationDriver,
    public G4ChordFinderDelegate<G4IntegrationDriver<G4BulirschStoer>>
{
  public:

    /**
     * Constructor for the concrete G4IntegrationDriver.
     *  @param[in] hminimum The minumum allowed step..
     *  @param[in] Boris Pointer to the Bulirsch-Stoer motion algorithm.
     *  @param[in] numberOfComponents The number of integration variables.
     *  @param[in] verbosity Flag for verbosity.
     */
    G4IntegrationDriver( G4double hminimum,
                         G4BulirschStoer* stepper,
                         G4int numberOfComponents = 6,
                         G4int statisticsVerbosity = 1);

    /**
     * Default Destructor.
     */
    ~G4IntegrationDriver() = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4IntegrationDriver(const G4IntegrationDriver&) = delete;
    G4IntegrationDriver& operator=(const G4IntegrationDriver&) = delete;

    /**
     * Computes the step to take, based on chord limits.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] chordDistance Maximum sagitta distance.
     *  @returns The length of step taken.
     */
    G4double AdvanceChordLimited(G4FieldTrack& track,
                                 G4double hstep,
                                 G4double eps,
                                 G4double chordDistance) override;

    /**
     * Dispatch interface method for initialisation/reset of driver.
     */
    void OnStartTracking() override;

    /**
     * Dispatch interface method for computing step. Does nothing here.
     */
    void OnComputeStep(const G4FieldTrack* track = nullptr) override;

    /**
     * The driver does not implement re-integration. Returns false.
     */
    G4bool DoesReIntegrate() const override;
   
    /**
     * Advances integration accurately by relative accuracy better than 'eps'.
     *  @param[in,out] track The current track in field.
     *  @param[in] stepLen Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] beginStep Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    G4bool AccurateAdvance( G4FieldTrack& track,
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
    G4bool QuickAdvance( G4FieldTrack& y_val,
                         const G4double dydx[],
                         G4double hstep,
                         G4double& missDist,
                         G4double& dyerr) override;

    /**
     * Takes one Step that is as large as possible while satisfying the
     * accuracy criterion.
     *  @param[in,out] y The current track state, y.
     *  @param[in] dydx dydx array.
     *  @param[in,out] curveLength Step start, x.
     *  @param[in] htry Step to attempt.
     *  @param[in] eps The relative accuracy.
     *  @param[out] hdid Step achieved.
     *  @param[out] hnext Proposed next step.
     */
    void OneGoodStep( G4double y[],
                      const G4double dydx[],
                      G4double& curveLength,
                      G4double htry,
                      G4double eps,
                      G4double& hdid,
                      G4double& hnext);

    /**
     * Getters for derivatives.
     */
    void GetDerivatives( const G4FieldTrack& track,
                               G4double dydx[]) const override;
    void GetDerivatives( const G4FieldTrack& track,
                               G4double dydx[],
                               G4double field[]) const override;

    /**
     * Setter and getter for verbosity.
     */
    void SetVerboseLevel(G4int level) override;
    G4int GetVerboseLevel() const override;

    /**
     * Computes the new step size .
     *  @param[in] errMaxNorm The normalised error.
     *  @param[in] hstepCurrent The current step size.
     *  @returns The new step size.
     */
    G4double ComputeNewStepSize(G4double errMaxNorm,
                                G4double hstepCurrent) override;

    /**
     * Getters and setter for the equation of motion.
     */
    G4EquationOfMotion* GetEquationOfMotion() override;
    const G4EquationOfMotion* GetEquationOfMotion() const;
    void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    /**
     * Getters for the stepper.
     */
    const G4MagIntegratorStepper* GetStepper() const override;
    G4MagIntegratorStepper* GetStepper() override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    void StreamInfo( std::ostream& os ) const override;
   
  private:

    G4int GetNumberOfVarialbles() const;

    G4double fMinimumStep;
    G4double fVerbosity;

    G4ModifiedMidpoint fMidpointMethod;
    G4BulirschStoer* bulirschStoer;

    G4double yIn[G4FieldTrack::ncompSVEC],
             yMid[G4FieldTrack::ncompSVEC],
             yMid2[G4FieldTrack::ncompSVEC],
             yOut[G4FieldTrack::ncompSVEC],
             yOut2[G4FieldTrack::ncompSVEC],
             yError[G4FieldTrack::ncompSVEC];


    G4double dydxCurrent[G4FieldTrack::ncompSVEC];
    G4double yCurrent[G4FieldTrack::ncompSVEC];

    G4double derivs[2][6][G4FieldTrack::ncompSVEC];

    const G4int interval_sequence[2];

    using ChordFinderDelegate =
          G4ChordFinderDelegate<G4IntegrationDriver<G4BulirschStoer>>;
};

#include "G4BulirschStoerDriver.icc"

#endif
