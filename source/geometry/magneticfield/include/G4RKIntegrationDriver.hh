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
// G4RKIntegrationDriver
//
// Class description:
//
// Driver class which controls the integration error of a 
// Runge-Kutta stepper 

// Created: D.Sorokin, 2017
// --------------------------------------------------------------------
#ifndef G4RKINTEGRATIONDRIVER_HH
#define G4RKINTEGRATIONDRIVER_HH

#include "G4VIntegrationDriver.hh"

template <class T>
class G4RKIntegrationDriver : public G4VIntegrationDriver
{
  public:

    G4RKIntegrationDriver(T* stepper);

    G4RKIntegrationDriver(const G4RKIntegrationDriver&) = delete;
    G4RKIntegrationDriver& operator=(const G4RKIntegrationDriver&) = delete;

    virtual void GetDerivatives(const G4FieldTrack& track,
                                      G4double dydx[]) const override;

    virtual void GetDerivatives(const G4FieldTrack& track,
                                      G4double dydx[],
                                      G4double field[]) const override;

    virtual G4double ComputeNewStepSize(G4double errMaxNorm, // normalised error
                                        G4double hstepCurrent) override final;
      // Taking the last step's normalised error, calculate
      // a step size for the next step.
      // - Limits the next step's size within a factor of the current one.

    virtual G4EquationOfMotion* GetEquationOfMotion() override;
    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual const T* GetStepper() const override;
    virtual T* GetStepper() override;

    virtual void  StreamInfo( std::ostream& os ) const override;
   
    // Accessors.
    G4double GetSafety() const;
    G4double GetPshrnk() const;
    G4double GetPgrow() const;

    virtual void RenewStepperAndAdjust(G4MagIntegratorStepper* stepper) override;

    void ReSetParameters(G4double safety = 0.9);
    void SetSafety(G4double valS);
      //  i) sets the exponents (pgrow & pshrnk),
      //     using the current Stepper's order,
      // ii) sets the safety

     G4int GetMaxNoSteps() const;
     void SetMaxNoSteps(G4int val);
       // Modify and Get the Maximum number of Steps that can be
       // taken for the integration of a single segment -
       // (ie a single call to AccurateAdvance).

     G4double GetSmallestFraction() const;
     void SetSmallestFraction(G4double val);

  protected:

    G4double ShrinkStepSize(G4double h, G4double error) const;
    G4double GrowStepSize(G4double h, G4double error) const;

    G4double ShrinkStepSize2(G4double h, G4double error2) const;
    G4double GrowStepSize2(G4double h, G4double error2) const;

    void UpdateErrorConstraints();

  private:

    inline void RenewStepperAndAdjustImpl(T* stepper);

    G4int fMaxNoSteps;

    G4int fMaxStepBase;
      // The (default) maximum number of steps is Base
      // divided by the order of Stepper

    // Parameters used to grow and shrink trial stepsize.
    //
    G4double safety;
    G4double pshrnk;   //  exponent for shrinking
    G4double pgrow;    //  exponent for growth

    // muximum error values for shrinking / growing (optimisation).
    //
    G4double errorConstraintShrink;
    G4double errorConstraintGrow;

    T* pIntStepper = nullptr;
};

#include "G4RKIntegrationDriver.icc"

#endif
