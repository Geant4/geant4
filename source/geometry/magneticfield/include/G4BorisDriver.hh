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
//   G4BorisDriver is a driver class using the second order Boris 
// method to integrate the equation of motion.
// 
//
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

    G4BorisDriver( G4double hminimum,
                   G4BorisScheme* Boris,
                   G4int numberOfComponents = 6,
                   bool verbosity = false);
   
    inline ~G4BorisDriver() = default;

    inline G4BorisDriver(const G4BorisDriver&) = delete;
    inline G4BorisDriver& operator=(const G4BorisDriver&) = delete;

    // 1. Core methods that advance the integration
    virtual G4bool AccurateAdvance( G4FieldTrack& track,
                                    G4double stepLen,
                                    G4double epsilon,
                                    G4double beginStep = 0) override;
    // Advance integration accurately - by relative accuracy better than 'epsilon'
   
     virtual G4bool QuickAdvance( G4FieldTrack& y_val,   // In/Out
                                  const G4double dydx[],    
                                  G4double       hstep,
                                  G4double&   missDist,  // Out: estimated sagitta
                                  G4double&   dyerr   ) override;
    // Attempt one integration step, and return estimated error 'dyerr'

    void OneGoodStep(G4double  yCurrentState[],  // In/Out: state ('y')
                     G4double& curveLength,      // In/Out: 'x'
                     G4double  htry,             // step to attempt
                     G4double  epsilon_rel,      // relative accuracy
                     G4double  restMass,         
                     G4double  charge,
                     G4double& hdid,             // Out: step achieved
                     G4double& hnext);           // Out: proposed next step
    //   Method to implement Accurate Advance
   
    // 2. Methods needed to co-work with G4ChordFinder
    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override
    {
      return ChordFinderDelegate::
             AdvanceChordLimitedImpl(track, hstep, eps, chordDistance);
    }

    virtual void OnStartTracking() override {
      ChordFinderDelegate::ResetStepEstimate();
    }

    virtual void OnComputeStep() override {};



    // 3. Does the method redo integrations when called to obtain values
    //      for internal, smaller intervals ? 
    //        (when needed to identify an intersection.) 
    virtual G4bool DoesReIntegrate() const override { return true; }
    //    It would be no if it just used interpolation to provide a result.

    // 4. Relevant for calculating a new step size to achieve required accuracy
    inline virtual G4double ComputeNewStepSize(
                          G4double  errMaxNorm,    // normalised error
                          G4double  hstepCurrent) override; // current step size

    G4double ShrinkStepSize2(G4double h, G4double error2) const;
    G4double GrowStepSize2(G4double h, G4double error2) const;
    // Calculate the next step size given the square of the relative error
   
    // 5. Auxiliary Methods ...
    virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[]) const override;

    virtual void GetDerivatives( const G4FieldTrack& track,
                                 G4double dydx[],
                                 G4double field[]) const override;

    inline virtual void SetVerboseLevel(G4int level) override;
    inline virtual G4int GetVerboseLevel() const override;

    inline virtual G4EquationOfMotion* GetEquationOfMotion() override;
    inline const G4EquationOfMotion* GetEquationOfMotion() const;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual void  StreamInfo( std::ostream& os ) const override;
     // Write out the parameters / state of the driver
   
    // 6. Not relevant for Boris and other non-RK methods
    inline virtual const G4MagIntegratorStepper* GetStepper() const override;
    inline virtual G4MagIntegratorStepper* GetStepper() override;

  private:
    inline G4int GetNumberOfVariables() const;

    inline void CheckStep(const G4ThreeVector& posIn,                              
                          const G4ThreeVector& posOut,
                          G4double hdid) const;
   
  private:
    // INVARIANTS -- remain unchanged during tracking / integration 
    // Parameters
    G4double fMinimumStep;
    bool     fVerbosity;

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
    static constexpr int      fMaxNoSteps = 300; 
    static constexpr G4double fSmallestFraction= 1e-12; // To avoid FP underflow !  ( 1.e-6 for single prec)

    static constexpr G4int    fIntegratorOrder= 2; //  2nd order method -- needed for error control
    static constexpr G4double fSafetyFactor = 0.9; //

    static constexpr G4double fMaxSteppingIncrease= 10.0; //  Increase no more than 10x   
    static constexpr G4double fMaxSteppingDecrease= 0.1;  //  Reduce   no more than 10x
    static constexpr G4double fPowerShrink = -1.0 / fIntegratorOrder;
    static constexpr G4double fPowerGrow   = -1.0 / (1.0 + fIntegratorOrder);

    static const G4double fErrorConstraintShrink;
    static const G4double fErrorConstraintGrow;
   
    using ChordFinderDelegate =
          G4ChordFinderDelegate<G4BorisDriver>;
};

#include "G4BorisDriver.icc"

#endif
