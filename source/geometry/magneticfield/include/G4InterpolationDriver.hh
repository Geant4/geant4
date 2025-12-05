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
// G4InterpolationDriver
//
// Class description:
//
// Driver class which uses Runge-Kutta stepper with interpolation property
// to integrate track with error control

// Author: Dmitry Sorokin (CERN), 26.09.2018
// --------------------------------------------------------------------
#ifndef G4INTERPOLATION_DRIVER_HH
#define G4INTERPOLATION_DRIVER_HH

#include "G4FieldUtils.hh"
#include "G4RKIntegrationDriver.hh"
#include "globals.hh"

#include <memory>
#include <vector>

/**
 * @brief G4InterpolationDriver is a templated driver class which uses
 * Runge-Kutta stepper with interpolation property to integrate track with
 * error control.
 */

template <class T, G4bool StepperCachesDchord = true>
class G4InterpolationDriver : public G4RKIntegrationDriver<T>
{
  public:

    /**
     * Constructor for G4IntegrationDriver.
     *  @param[in] hminimum Minimum allowed step.
     *  @param[in] stepper Pointer to the stepper algorithm.
     *  @param[in] numberOfComponents The number of integration variables,
     *             if not matching stepper's number of variables, issue exception.
     *  @param[in] statisticsVerbosity Verbosity level.
     */
    G4InterpolationDriver(G4double hminimum,
                          T* stepper,
                          G4int numberOfComponents = 6,
                          G4int statisticsVerbosity = 0);

    /**
     * Destructor. Provides statistics if verbosity level is greater than zero.
     */
    ~G4InterpolationDriver() override;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4InterpolationDriver(const G4InterpolationDriver&) = delete;
    const G4InterpolationDriver& operator=(const G4InterpolationDriver&) = delete;

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
    void OnComputeStep(const G4FieldTrack* /*track*/ = nullptr) override;

    /**
     * The driver does not implement re-integration. Returns false.
     */
    G4bool DoesReIntegrate() const override { return false; }

    /**
     * Advances integration accurately by relative accuracy better than 'eps'.
     * On output the track is replaced by the value at the end of interval.
     *  @param[in,out] track The current track in field.
     *  @param[in] hstep Proposed step length.
     *  @param[in] eps Requested accuracy, y_err/hstep.
     *  @param[in] hinitial Initial minimum integration step.
     *  @returns true if integration succeeds.
     */
    G4bool AccurateAdvance(G4FieldTrack& track,
                           G4double hstep,
                           G4double eps,  // Requested y_err/hstep
                           G4double hinitial = 0) override;

    /**
     * Setter and getter for verbosity.
     */
    void SetVerboseLevel(G4int level) override;
    G4int GetVerboseLevel() const override;

    /**
     * Writes out to stream the parameters/state of the driver.
     */
    void StreamInfo(std::ostream& os) const override;

  protected:

    struct InterpStepper
    {
      std::unique_ptr<T> stepper;
      G4double begin;
      G4double end;
      G4double inverseLength;
    };

    using StepperIterator = typename std::vector<InterpStepper>::iterator;
    using ConstStepperIterator = typename std::vector<InterpStepper>::const_iterator;

    /**
     * Takes one Step that is as large as possible while satisfying the
     * accuracy criterion.
     *  @param[in] it Stepper iterator.
     *  @param[in,out] y The current track state, y.
     *  @param[in] dydx dydx array.
     *  @param[in,out] hstep Step to attempt.
     *  @param[in] eps The relative accuracy.
     *  @param[in] curveLength Step start, x.
     *  @param[in,out] track Pointer to the Field track. Not used.
     *  @returns The step achieved.
     */
    virtual G4double OneGoodStep(StepperIterator it,
                                 field_utils::State& y,
                                 field_utils::State& dydx,
                                 G4double& hstep,
                                 G4double eps,
                                 G4double curveLength,
                                 G4FieldTrack* track = nullptr);

    /**
     * Track interpolation.
     *  @param[in] curveLength Step start, x.
     *  @param[in,out] y The current track state, y.
     */
    void Interpolate(G4double curveLength, field_utils::State& y) const;

    /**
     * Wrapper method for interpolation.
     */
    void InterpolateImpl(G4double curveLength,
                         ConstStepperIterator it,
                         field_utils::State& y) const;

    /**
     * Methods for calculation of chord step and distance.
     */
    G4double DistChord(const field_utils::State& yBegin,
                             G4double curveLengthBegin,
                       const field_utils::State& yEnd,
                             G4double curveLengthEnd) const;
    G4double FindNextChord(const field_utils::State& yBegin,
                                 G4double curveLengthBegin,
                                 field_utils::State& yEnd,
                                 G4double curveLengthEnd,
                                 G4double dChord,
                                 G4double maxChordDistance);
    G4double CalcChordStep(G4double stepTrialOld,
                           G4double dChordStep,
                           G4double fDeltaChord);


    /**
     * Internal methods for printing/checking the state.
     */
    void PrintState() const;
    void CheckState() const;

    /**
     * Increments number of trials and calls.
     */
    void AccumulateStatistics(G4int noTrials);

  protected:

    std::vector<InterpStepper> fSteppers;
    StepperIterator fLastStepper;
    G4bool fKeepLastStepper = false;

    /** Memory of last good step size for integration. */
    G4double fhnext = DBL_MAX;

    /** Minimum Step allowed (in units of length). */
    G4double fMinimumStep;

    G4double fChordStepEstimate = DBL_MAX;
    const G4double fFractionNextEstimate = 0.98;  // Constant
    const G4double fSmallestCurveFraction = 0.01;  // Constant

    G4int fVerboseLevel;  // Parameter

    field_utils::State fdydx;
    G4bool fFirstStep = true;

    G4int fMaxTrials = 100;  // Constant
    G4int fTotalStepsForTrack = 0;

    /** Statistics. */
    G4int fTotalNoTrials = 0;
    G4int fNoCalls = 0;
    G4int fmaxTrials = 0;

    using Base = G4RKIntegrationDriver<T>;
};

#include "G4InterpolationDriver.icc"

#endif
