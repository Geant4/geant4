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
// G4DormandPrince745
//
// Class desription:
//
//  An implementation of the 5th order embedded RK method from the paper:
//  J. R. Dormand and P. J. Prince, "A family of embedded Runge-Kutta formulae"
//  Journal of computational and applied Math., vol.6, no.1, pp.19-26, 1980.
//
//  DormandPrince7 - 5(4) embedded RK method
//

// Created: Somnath Banerjee, Google Summer of Code 2015, 25 May 2015
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4DORMAND_PRINCE_745_HH
#define G4DORMAND_PRINCE_745_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldUtils.hh"

class G4DormandPrince745 : public G4MagIntegratorStepper
{
  public:

    G4DormandPrince745(G4EquationOfMotion* equation,
                       G4int numberOfVariables = 6);

    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[]) override;

    void Stepper(const G4double yInput[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOutput[],
                       G4double yError[],
                       G4double dydxOutput[]);

    inline void SetupInterpolation() {}

    inline void Interpolate(G4double tau, G4double yOut[]) const
    {
      Interpolate4thOrder(yOut, tau);
    }
      // For calculating the output at the tau fraction of Step

    G4double DistChord() const override;

    G4int IntegratorOrder() const override { return 4; }

    const G4String& StepperType() const        { return gStepperType; }
    const G4String& StepperDescription() const { return gStepperDescription; }
   
    const field_utils::State& GetYOut() const { return fyOut; }

    void Interpolate4thOrder(G4double yOut[], G4double tau) const;

    void SetupInterpolation5thOrder();
    void Interpolate5thOrder(G4double yOut[], G4double tau) const;

    G4EquationOfMotion* GetSpecificEquation() { return GetEquationOfMotion(); }

  private:

    static const G4String gStepperType;
    static const G4String gStepperDescription;
      // Name and description of this steppers
      // plus details of its implementation

    field_utils::State ak2, ak3, ak4, ak5, ak6, ak7, ak8, ak9;
    field_utils::State fyIn, fyOut, fdydxIn;

    G4double fLastStepLength = -1.0;
};

#endif
