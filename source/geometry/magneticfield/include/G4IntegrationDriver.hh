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
// $Id: G4IntegrationDriver.hh 106739 2017-10-20 14:45:37Z dmsoroki $
//
//
// class G4IntegrationDriver
//
// Class description:
//
// Driver class which controls the integration error of a 
// Runge-Kutta stepper 

// History:
// - Created. D.Sorokin
// --------------------------------------------------------------------

#ifndef G4IntegrationDriver_HH
#define G4IntegrationDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4VIntegrationDriver.hh"

template <class T>
class G4IntegrationDriver : public G4VIntegrationDriver {
public:
    G4IntegrationDriver(  G4double hminimum,
                          T*       stepper,
                          G4int    numberOfComponents = 6,
                          G4int    statisticsVerbosity = 1);

    virtual ~G4IntegrationDriver() override;

    G4IntegrationDriver(const G4IntegrationDriver &) = delete;
    const G4IntegrationDriver& operator =(const G4IntegrationDriver &) = delete;

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
        G4double& dchord_step,
        G4double& dyerr) override;

    virtual void GetDerivatives(
         const G4FieldTrack &track,
         G4double dydx[]) const override;

    // Taking the last step's normalised error, calculate
    // a step size for the next step.
    // Do not limit the next step's size within a factor of the
    // current one.
    virtual G4double ComputeNewStepSize(
        G4double errMaxNorm,             // normalised error
        G4double hstepCurrent) override; // current step size

    virtual void SetVerboseLevel(G4int newLevel) override;
    virtual G4int GetVerboseLevel() const override;

    virtual G4EquationOfMotion* GetEquationOfMotion() override;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual const T* GetStepper() const override;
    virtual T* GetStepper() override;

    // Accessors.
    G4double GetMinimumStep() const;
    G4double GetSafety() const;
    G4double GetPshrnk() const;
    G4double GetPgrow() const;

    // Sets a new stepper pItsStepper for this driver. Then it calls
    // ReSetParameters to reset its parameters accordingly.
    void RenewStepperAndAdjust(T *pItsStepper);

    //  i) sets the exponents (pgrow & pshrnk),
    //     using the current Stepper's order,
    // ii) sets the safety
    void ReSetParameters(G4double safety = 0.9);

    void SetMinimumStep(G4double newval);
    void SetSafety(G4double valS);

    // This takes one Step that is of size htry, or as large 
    // as possible while satisfying the accuracy criterion of:
    //     yerr < eps * |y_end-y_start|
    void OneGoodStep(      G4double  yVar[],  // InOut
                     const G4double  dydx[],
                           G4double& curveLength,
                           G4double  htry,
                           G4double  eps,
                           G4double& hdid,
                           G4double& hnext);

     // Modify and Get the Maximum number of Steps that can be
     // taken for the integration of a single segment -
     // (ie a single call to AccurateAdvance).
     G4int GetMaxNoSteps() const;
     void SetMaxNoSteps( G4int val);

     G4double GetSmallestFraction() const;
     void SetSmallestFraction(G4double val);

private:
     G4double ShrinkStepSize(G4double h, G4double error) const;
     G4double GrowStepSize(G4double h, G4double error) const;
     void UpdateErrorConstraints();

     void CheckStep(
         const G4ThreeVector& posIn, const G4ThreeVector& posOut, G4double hdid);

     // Minimum Step allowed in a Step (in absolute units)
     G4double  fMinimumStep;

     // Smallest fraction of (existing) curve length - in relative units
     // below this fraction the current step will be the last
     G4double  fSmallestFraction;
     //  Expected range: smaller than 0.1 * epsilon and bigger than 5e-13
     //    ( Note: this range is not enforced. )
     G4int   fMaxNoSteps;

     // The (default) maximum number of steps is Base
     // divided by the order of Stepper
     static constexpr G4int  fMaxStepBase = 250;

     // Parameters used to grow and shrink trial stepsize.
     G4double safety;
     G4double pshrnk;   //  exponent for shrinking
     G4double pgrow;    //  exponent for growth

     // muximum error values for shrinking / growing (optimisation).
     G4double errorConstraintShrink;
     G4double errorConstraintGrow;

     // Maximum stepsize increase/decrease factors.
     static constexpr G4double max_stepping_increase = 5;
     static constexpr G4double max_stepping_decrease = 0.1;

     T* pIntStepper;

     // Step Statistics
     unsigned long fNoTotalSteps, fNoBadSteps, fNoGoodSteps;

     G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues

     G4int fNoQuickAvanceCalls;
};

#include "G4IntegrationDriver.icc"

#endif
