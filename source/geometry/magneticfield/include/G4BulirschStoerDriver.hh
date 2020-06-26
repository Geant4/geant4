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
// G4IntegrationDriver<G4BulirschStoer> is a driver class using
// Bulirsch-Stoer method to integrate the equation of motion.

// Author: Dmitry Sorokin, Google Summer of Code 2016
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4BULIRSCH_STOER_DRIVER_HH
#define G4BULIRSCH_STOER_DRIVER_HH

#include "G4IntegrationDriver.hh"
#include "G4BulirschStoer.hh"
#include "G4ChordFinderDelegate.hh"

template <>
class G4IntegrationDriver<G4BulirschStoer>: 
    public G4VIntegrationDriver,
    public G4ChordFinderDelegate<G4IntegrationDriver<G4BulirschStoer>>
{
  public:

    G4IntegrationDriver( G4double hminimum,
                         G4BulirschStoer* stepper,
                         G4int numberOfComponents = 6,
                         G4int statisticsVerbosity = 1);

    ~G4IntegrationDriver() = default;

    G4IntegrationDriver(const G4IntegrationDriver&) = delete;
    G4IntegrationDriver& operator=(const G4IntegrationDriver&) = delete;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override
    {
      return ChordFinderDelegate::
             AdvanceChordLimitedImpl(track, hstep, eps, chordDistance);
    }

    virtual void OnStartTracking() override
    {
      ChordFinderDelegate::ResetStepEstimate();
    }

    virtual void OnComputeStep() override {};

    virtual G4bool DoesReIntegrate() const override { return false; }  /// ????
   
    virtual G4bool AccurateAdvance( G4FieldTrack& track,
                                    G4double stepLen,
                                    G4double eps,
                                    G4double beginStep = 0) override;

    virtual G4bool QuickAdvance( G4FieldTrack& y_val,
                                 const G4double dydx[],
                                 G4double hstep,
                                 G4double& missDist,
                                 G4double& dyerr) override;

    void OneGoodStep( G4double y[],
                      const G4double dydx[],
                      G4double& curveLength,
                      G4double htry,
                      G4double eps,
                      G4double& hdid,
                      G4double& hnext);

    virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[]) const override;

    virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[],
                                 G4double field[]) const override;

    virtual void SetVerboseLevel(G4int level) override;
    virtual G4int GetVerboseLevel() const override;

    virtual G4double ComputeNewStepSize(
                          G4double  errMaxNorm,    // normalised error
                          G4double  hstepCurrent) override; // current step size

    virtual G4EquationOfMotion* GetEquationOfMotion() override;
    const G4EquationOfMotion* GetEquationOfMotion() const;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual const G4MagIntegratorStepper* GetStepper() const override;
    virtual G4MagIntegratorStepper* GetStepper() override;

    virtual void  StreamInfo( std::ostream& os ) const override;
     // Write out the parameters / state of the driver
   
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
