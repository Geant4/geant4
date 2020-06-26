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

// Created: D.Sorokin, 2018
// --------------------------------------------------------------------
#ifndef G4INTERPOLATION_DRIVER_HH
#define G4INTERPOLATION_DRIVER_HH

#include "G4RKIntegrationDriver.hh"
#include "G4FieldUtils.hh"

#include "globals.hh"

#include <vector>
#include <memory>

template <class T>
class G4InterpolationDriver : public G4RKIntegrationDriver<T>
{
  public:

    G4InterpolationDriver(G4double hminimum,
                          T* stepper,
                          G4int numberOfComponents = 6,
                          G4int statisticsVerbosity = 0);

    virtual ~G4InterpolationDriver() override;

    G4InterpolationDriver(const G4InterpolationDriver&)= delete;
    const G4InterpolationDriver& operator=(const G4InterpolationDriver&)= delete;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override;

    virtual void OnStartTracking() override;
    virtual void OnComputeStep() override;
    virtual G4bool DoesReIntegrate() const override { return false; }
     // Interpolation driver does not recalculate when AccurateAdvance is called
     //  -- reintegration would require other calls
   
    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps,    // Requested y_err/hstep
                                   G4double hinitial = 0) override;
      // Integrates ODE from current s (s=s0) to s=s0+h with accuracy eps.
      // On output track is replaced by value at end of interval.
      // The concept is similar to the odeint routine from NRC p.721-722.

    virtual void SetVerboseLevel(G4int level) override;
    virtual G4int GetVerboseLevel() const override;

    virtual void  StreamInfo( std::ostream& os ) const override;
   
  private:

    struct InterpStepper
    {
        std::unique_ptr<T> stepper;
        G4double begin;
        G4double end;
        G4double inverseLength;
    };

    using StepperIterator = typename std::vector<InterpStepper>::iterator;
    using ConstStepperIterator = typename std::vector<InterpStepper>::const_iterator;

    G4double OneGoodStep(StepperIterator it,
                         field_utils::State& y,
                         field_utils::State& dydx,
                         G4double& hstep,
                         G4double eps,
                         G4double curveLength);
      // This takes one Step that is of size htry, or as large 
      // as possible while satisfying the accuracy criterion of:
      //     yerr < eps * |y_end-y_start|
      // return hdid

    void Interpolate(G4double curveLength, field_utils::State& y) const;

    void InterpolateImpl(G4double curveLength, 
                         ConstStepperIterator it, 
                         field_utils::State& y) const;

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

    void PrintState() const;

    void CheckState() const;

    void AccumulateStatistics(G4int noTrials);

  private:

    std::vector<InterpStepper> fSteppers;
    StepperIterator fLastStepper;
    G4bool fKeepLastStepper = false;

    G4double fhnext = DBL_MAX;  // Memory of last good step size for integration

    // Minimum Step allowed in a Step (in units of length)   // Parameter
    G4double fMinimumStep;

    G4double fChordStepEstimate = DBL_MAX;
    const G4double fFractionNextEstimate = 0.98;             // Constant
    const G4double fSmallestCurveFraction = 0.01;            // Constant

    G4int fVerboseLevel;                                     // Parameter

    field_utils::State fdydx;
    G4bool fFirstStep = true;

    const G4int fMaxTrials = 100;                            // Constant
    G4int fTotalStepsForTrack = 0;

    // statistics
    G4int fTotalNoTrials = 0;
    G4int fNoCalls = 0;
    G4int fmaxTrials = 0;

    using Base = G4RKIntegrationDriver<T>;
};

#include "G4InterpolationDriver.icc"

#endif
