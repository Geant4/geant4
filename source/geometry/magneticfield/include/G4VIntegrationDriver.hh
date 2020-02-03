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
// G4VIntegrationDriver
//
// Class description:
//
// Abstract base class for 'driver' classes which are responsible for
// undertaking integration of an state given an equation of motion and
// within acceptable error bound(s).
//
// Different integration methods are meant to be provided via this
// common interface, and can span the original type (explicit Runge Kutta
// methods), enhanced RK methods and alternatives such as the
// Bulirsch-Stoer and multi-step methods.
//
// The drivers' key mission is to insure that the error is below set values. 

// Author: Dmitry Sorokin, Google Summer of Code 2017
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4VINTEGRATION_DRIVER_HH
#define G4VINTEGRATION_DRIVER_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4MagIntegratorStepper.hh"

class G4VIntegrationDriver
{
  public:

    virtual ~G4VIntegrationDriver() = default;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) = 0;

    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps, // Requested y_err/hstep
                                   G4double hinitial = 0) = 0;

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) = 0;
    virtual G4EquationOfMotion* GetEquationOfMotion() = 0;

    virtual void RenewStepperAndAdjust(G4MagIntegratorStepper* pItsStepper);
      // Method for compatibility -- relevant only for G4MagIntegratorDriver
   
    virtual void SetVerboseLevel(G4int level) = 0;
    virtual G4int GetVerboseLevel() const = 0;

    virtual void OnComputeStep() = 0;

    virtual void OnStartTracking() = 0;

  public:  // without description

    //[[deprecated("will be removed")]]
    virtual G4bool QuickAdvance(G4FieldTrack& /*track*/,   // INOUT
                                const G4double /*dydx*/[],
                                G4double /*hstep*/,
                                G4double& /*dchord_step*/,
                                G4double& /*dyerr*/) { return false; }

    //[[deprecated("will be removed")]]
    virtual void GetDerivatives(const G4FieldTrack& track,
                                G4double dydx[]) const = 0;

    //[[deprecated("will be removed")]]
    virtual void GetDerivatives(const G4FieldTrack& track,
                                G4double dydx[],
                                G4double field[]) const = 0;

    //[[deprecated("use GetEquationOfMotion() instead of GetStepper()->GetEquationOfMotion()")]]
    virtual const G4MagIntegratorStepper* GetStepper() const = 0;
    virtual G4MagIntegratorStepper* GetStepper() = 0;

    //[[deprecated("will be removed")]]
    virtual G4double ComputeNewStepSize(G4double errMaxNorm, // normalised error
                                        G4double hstepCurrent) = 0;
      // Taking the last step's normalised error, calculate
      // a step size for the next step.
      // Do not limit the next step's size within a factor of the current one.

    virtual G4bool DoesReIntegrate() = 0;
      // Whether the driver implementates re-integration
      //  - original Integration driver will re-start and re-calculate interval => yes
      //  - Interpolation Driver does not recalculate (it interpolates)
      // Basically answer: does this driver *Recalculate* when AccurateAdvance is called ?
  protected:

    static constexpr G4double max_stepping_increase = 5;
    static constexpr G4double max_stepping_decrease = 0.1;
};

#endif
