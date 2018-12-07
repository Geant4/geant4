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
// class G4FSALIntegrationDriver
//
// Class description:
//
// Driver class which controls the integration error of a 
// Runge-Kutta stepper 

// History:
// - Created. D.Sorokin
// --------------------------------------------------------------------

#ifndef G4FSALIntegrationDriver_HH
#define G4FSALIntegrationDriver_HH

#include "G4RKIntegrationDriver.hh"
#include "G4ChordFinderDelegate.hh"

template <class T>
class G4FSALIntegrationDriver : public G4RKIntegrationDriver<T> ,
                                public G4ChordFinderDelegate<G4FSALIntegrationDriver<T>> {
public:
    G4FSALIntegrationDriver(  G4double hminimum,
                          T*       stepper,
                          G4int    numberOfComponents = 6,
                          G4int    statisticsVerbosity = 1);

    virtual ~G4FSALIntegrationDriver() override;

    G4FSALIntegrationDriver(const G4FSALIntegrationDriver &) = delete;
    const G4FSALIntegrationDriver& operator =(const G4FSALIntegrationDriver &) = delete;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override
    {
        return ChordFinderDelegate::AdvanceChordLimitedImpl(track, hstep, eps, chordDistance);
    }

    virtual void OnStartTracking() override
    {
        ChordFinderDelegate::ResetStepEstimate();
    }

    virtual void OnComputeStep() override {};

    // Integrates ODE from current s (s=s0) to s=s0+h with accuracy eps.
    // On output track is replaced by value at end of interval.
    // The concept is similar to the odeint routine from NRC p.721-722.
    virtual G4bool AccurateAdvance(
        G4FieldTrack& track,
        G4double hstep,
        G4double eps,                     // Requested y_err/hstep
        G4double hinitial = 0) override;  // Suggested 1st interval

    // QuickAdvance just tries one Step - it does not ensure accuracy.
    virtual G4bool QuickAdvance(
        G4FieldTrack& fieldTrack,
        const G4double dydx[],
        G4double hstep,
        G4double inverseCurvatureRadius,
        G4double& dchord_step,
        G4double& dyerr) override;

    virtual void SetVerboseLevel(G4int newLevel) override;
    virtual G4int GetVerboseLevel() const override;

    // Accessors.
    G4double GetMinimumStep() const;
    void SetMinimumStep(G4double newval);

    // This takes one Step that is of size htry, or as large 
    // as possible while satisfying the accuracy criterion of:
    //     yerr < eps * |y_end-y_start|
    void OneGoodStep(G4double y[],  // InOut
                     G4double dydx[],
                     G4double& curveLength,
                     G4double htry,
                     G4double eps,
                     G4double& hdid,
                     G4double& hnext);

     G4double GetSmallestFraction() const;
     void SetSmallestFraction(G4double val);

protected:
    void IncrementQuickAdvanceCalls();

private:
    void CheckStep(const G4ThreeVector& posIn, 
                   const G4ThreeVector& posOut, 
                   G4double hdid);

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

    G4int fNoQuickAvanceCalls;
    G4int fNoAccurateAdvanceCalls;
    G4int fNoAccurateAdvanceBadSteps;
    G4int fNoAccurateAdvanceGoodSteps;

    using Base = G4RKIntegrationDriver<T>;
    using ChordFinderDelegate = G4ChordFinderDelegate<G4FSALIntegrationDriver<T>>;
};

#include "G4FSALIntegrationDriver.icc"

#endif
