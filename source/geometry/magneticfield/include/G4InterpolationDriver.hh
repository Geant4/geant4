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
//
//
//
// class G4InterpolationDriver
//
// Class description:
//
// Driver class which uses Runge-Kutta stepper with interpolation property
// to integrate track with error control 

// History:
// - Created. D.Sorokin
// --------------------------------------------------------------------

#ifndef G4InterpolationDriver_HH
#define G4InterpolationDriver_HH

#include "G4RKIntegrationDriver.hh"

using State = G4double[G4FieldTrack::ncompSVEC];

template <class T>
class G4InterpolationDriver : public G4RKIntegrationDriver<T> {
public:
    G4InterpolationDriver(  G4double hminimum,
                          T*       stepper,
                          G4int    numberOfComponents = 6,
                          G4int    statisticsVerbosity = 1);

    virtual ~G4InterpolationDriver() override;

    G4InterpolationDriver(const G4InterpolationDriver &) = delete;
    const G4InterpolationDriver& operator =(const G4InterpolationDriver &) = delete;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override;

    virtual void OnStartTracking() override
    {
        fhnext = 0;
        fChordStepEstimate = 0;
    };

    // Integrates ODE from current s (s=s0) to s=s0+h with accuracy eps.
    // On output track is replaced by value at end of interval.
    // The concept is similar to the odeint routine from NRC p.721-722.
    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps,                     // Requested y_err/hstep
                                   G4double hinitial = 0) override;  // Suggested 1st interval

    virtual void SetVerboseLevel(G4int newLevel) override;
    virtual G4int GetVerboseLevel() const override;

    virtual void OnComputeStep() override 
    { 
        fIntegrationInterval = {DBL_MAX, -DBL_MAX}; 
    }

    // Accessors.
    G4double GetMinimumStep() const;
    void SetMinimumStep(G4double newval);

    // This takes one Step that is of size htry, or as large 
    // as possible while satisfying the accuracy criterion of:
    //     yerr < eps * |y_end-y_start|
    void OneGoodStep(const G4double  yVar[],  // InOut
                     const G4double  dydx[],
                           G4double  htry,
                           G4double  eps,
                           G4double& hdid,
                           G4double& hnext);

     G4double GetSmallestFraction() const;
     void SetSmallestFraction(G4double val);

private:
    G4double FindNextChord(State& y,
                           G4double hstart,
                           G4double hmax,
                           G4double chordDistance);
    G4double BinsearchChord(State& y,
                           G4double hstart,
                           G4double hmax,
                           G4double chordDistance);


    void CheckStep(const G4ThreeVector& posIn, 
                   const G4ThreeVector& posOut, 
                   G4double hdid);

    std::pair<G4double, G4double> fIntegrationInterval;
    G4double fhnext;

    // Minimum Step allowed in a Step (in absolute units)
    G4double fMinimumStep;

    // Smallest fraction of (existing) curve length - in relative units
    // below this fraction the current step will be the last
    G4double fSmallestFraction;
    //  Expected range: smaller than 0.1 * epsilon and bigger than 5e-13
    //    ( Note: this range is not enforced. )

    // Verbosity level for printing (debug, ..)
    // Could be varied during tracking - to help identify issues
    G4int fVerboseLevel;

    G4int fNoAdvanceChordLimitedCalls;
    G4int fNoAdvanceChordLimitedSmallSteps;
    G4int fNoAdvanceChordLimitedFullSteps;
    G4int fNoAccurateAdvanceCalls;    
    G4int fNoAccurateAdvanceBadSteps;
    G4int fNoAccurateAdvanceGoodSteps;
    G4int fMaxTrials;

    G4double fChordStepEstimate = 0;

    using Base = G4RKIntegrationDriver<T>;
};

#include "G4InterpolationDriver.icc"

#endif
