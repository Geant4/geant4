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
// G4Boris
//
// Class description:
//
// G4BorisDriver is a driver class using
// boris method to integrate the equation of motion.

// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun 
// --------------------------------------------------------------------
#ifndef G4BORIS_DRIVER_HH
#define G4BORIS_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4BorisScheme.hh"
#include "G4ChordFinderDelegate.hh"


class G4BorisDriver: 
    public G4VIntegrationDriver,
    public G4ChordFinderDelegate<G4BorisDriver>
{
  public:

    inline G4BorisDriver( G4double hminimum,
                         G4BorisScheme* Boris,
                         G4int numberOfComponents = 6,
                         G4int statisticsVerbosity = 1);

    inline ~G4BorisDriver() = default;

    inline G4BorisDriver(const G4BorisDriver&) = delete;
   inline  G4BorisDriver& operator=(const G4BorisDriver&) = delete;

    inline virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override
    {
      return ChordFinderDelegate::
             AdvanceChordLimitedImpl(track, hstep, eps, chordDistance);
    }

    inline virtual void OnStartTracking() override
    {
      ChordFinderDelegate::ResetStepEstimate();
    }

    inline virtual void OnComputeStep() override {};

   inline  virtual G4bool DoesReIntegrate() const override { return false; }  /// ????
   
    inline virtual G4bool AccurateAdvance( G4FieldTrack& track,
                                    G4double stepLen,
                                    G4double eps,
                                    G4double beginStep = 0) override;

    inline virtual G4bool QuickAdvance( G4FieldTrack& y_val,
                                 const G4double dydx[],
                                 G4double hstep,
                                 G4double& missDist,
                                 G4double& dyerr) override;

    // inline void OneGoodStep(G4FieldTrack& track,
    //                    G4double y[],
    //                   const G4double dydx[],
    //                   G4double& curveLength,
    //                   G4double htry,
    //                   G4double eps,
    //                   G4double& hdid,
    //                   G4double& hnext);

    inline virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[]) const override;

    inline virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[],
                                 G4double field[]) const override;

    inline virtual void SetVerboseLevel(G4int level) override;
    inline virtual G4int GetVerboseLevel() const override;

    inline virtual G4double ComputeNewStepSize(
                          G4double  errMaxNorm,    // normalised error
                          G4double  hstepCurrent) override; // current step size

    inline virtual G4EquationOfMotion* GetEquationOfMotion() override;
    inline const G4EquationOfMotion* GetEquationOfMotion() const;
    inline virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

   inline  virtual const G4MagIntegratorStepper* GetStepper() const override;
    inline virtual G4MagIntegratorStepper* GetStepper() override;

    inline virtual void  StreamInfo( std::ostream& os ) const override;
     // Write out the parameters / state of the driver
   
  private:

    inline G4int GetNumberOfVarialbles() const;

    G4double fMinimumStep;
    G4double fVerbosity;

    G4BorisScheme* boris;

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
          G4ChordFinderDelegate<G4BorisDriver>;
};

#include "G4BorisDriver.icc"

#endif
