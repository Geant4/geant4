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
// G4FSALIntegrationDriver
//
// Class description:
//
// Driver class which controls the integration error of a Runge-Kutta stepper 

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2017), 20.10.2017
// --------------------------------------------------------------------
#ifndef G4FSALINTEGRATIONDRIVER_HH
#define G4FSALINTEGRATIONDRIVER_HH

#include "G4RKIntegrationDriver.hh"
#include "G4ChordFinderDelegate.hh"

/**
 * @brief G4FSALIntegrationDriver is a templated driver class which controls
 * the integration error of a Runge-Kutta stepper.
 */

template <class T>
class G4FSALIntegrationDriver : public G4RKIntegrationDriver<T>,
                                public G4ChordFinderDelegate<G4FSALIntegrationDriver<T>>
{
  public:

    /**
     * Constructor for G4FSALIntegrationDriver.
     *  @param[in] hminimum Minimum allowed step.
     *  @param[in] stepper Pointer to the stepper algorithm.
     *  @param[in] numberOfComponents The number of integration variables,
     *             if not matching stepper's number of variables, issue exception.
     *  @param[in] statisticsVerbosity Verbosity level.
     */
    inline G4FSALIntegrationDriver(G4double hminimum,
                                   T* stepper,
                                   G4int numberOfComponents = 6,
                                   G4int statisticsVerbosity = 1);

    /**
     * Destructor. Provides statistics if verbosity level is greater than zero.
     */
    inline ~G4FSALIntegrationDriver() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4FSALIntegrationDriver(const G4FSALIntegrationDriver&) = delete;
    G4FSALIntegrationDriver& operator=(const G4FSALIntegrationDriver&) = delete;

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
    inline void OnComputeStep(const G4FieldTrack* /*track*/ = nullptr) override;

    /**
     * The driver does implement re-integration. Returns true.
     */
    inline G4bool DoesReIntegrate() const override;
   
    /**
     * Advances integration accurately by relative accuracy better than 'eps'.
     * On output the track is replaced by the value at the end of interval.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] hinitial Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    inline G4bool AccurateAdvance(G4FieldTrack& track,
                                  G4double hstep,
                                  G4double eps, // Requested y_err/hstep
                                  G4double hinitial = 0.0) override;

    /**
     * Attempts one integration step, and returns estimated error 'dyerr'.
     * It does not ensure accuracy.
     *  @param[in,out] fieldTrack The current track in field.
     *  @param[in] dydx dydx array.
     *  @param[in] hstep Proposed step length.
     *  @param[out] dchord_step Estimated sagitta distance.
     *  @param[out] dyerr Estimated error.
     *  @returns true if integration succeeds.
     */
    inline G4bool QuickAdvance(G4FieldTrack& fieldTrack,
                               const G4double dydx[],
                               G4double hstep,
                               G4double& dchord_step,
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
    inline void OneGoodStep(G4double y[],  // InOut
                            G4double dydx[],
                            G4double& curveLength,
                            G4double htry,
                            G4double eps,
                            G4double& hdid,
                            G4double& hnext);

    /**
     * Setter and getter for verbosity.
     */
    inline void SetVerboseLevel(G4int newLevel) override;
    inline G4int GetVerboseLevel() const override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    inline void StreamInfo( std::ostream& os ) const override;
   
    /**
     * Getter and Setter for minimum allowed step.
     */
    inline G4double GetMinimumStep() const;
    inline void SetMinimumStep(G4double newval);

    /**
     * Getter and Setter for smallest fraction.
     */
    inline G4double GetSmallestFraction() const;
    inline void SetSmallestFraction(G4double val);

  protected:

    /**
     * Increments the counter for the number of calls to QuickAdvance().
     */
    inline void IncrementQuickAdvanceCalls();

  private:

    /**
     * Checks accuracy of step distance on the end point.
     */
    inline void CheckStep(const G4ThreeVector& posIn, 
                          const G4ThreeVector& posOut, G4double hdid);

  private:

    /** Minimum Step allowed in a Step (in absolute units). */
    G4double fMinimumStep;

    /** Smallest fraction of (existing) curve length in relative units.
     *  Below this fraction the current step will be the last.
     *  The expected range: smaller than 0.1 * epsilon and bigger than 5e-13
     *  (range not enforced). */
    G4double fSmallestFraction{1e-12};

    /** Verbosity level for printing (debug, etc..)
     *  Could be varied during tracking to help identifying issues. */
    G4int fVerboseLevel;

    G4int fNoQuickAvanceCalls{0};
    G4int fNoAccurateAdvanceCalls{0};
    G4int fNoAccurateAdvanceBadSteps{0};
    G4int fNoAccurateAdvanceGoodSteps{0};

    using Base = G4RKIntegrationDriver<T>;
    using ChordFinderDelegate = G4ChordFinderDelegate<G4FSALIntegrationDriver<T>>;
};

#include "G4FSALIntegrationDriver.icc"

#endif
