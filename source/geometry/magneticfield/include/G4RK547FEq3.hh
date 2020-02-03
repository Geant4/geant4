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
// G4RK547FEq3
//
// Class description:
//
// An implementation of the 7 stage embedded Runge-Kutta 4,5 pair (RK547FEq3)
// from the paper:
//    D. J. Higham and G. Hall,
//    "Embedded Runge-Kutta formulae with stable equilibrium states",
//    J. Comput. Appl. Math., vol. 29, no. 1, pp. 25-33, 1990.

// Author: Dmitry Sorokin, Google Summer of Code 2017
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4RK547FEq3_HH
#define G4RK547FEq3_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

class G4RK547FEq3 : public G4MagIntegratorStepper
{
  public:

    G4RK547FEq3(G4EquationOfMotion* EqRhs, G4int integrationVariables = 6);

    G4RK547FEq3(const G4RK547FEq3&) = delete;
    G4RK547FEq3& operator = (const G4RK547FEq3&) = delete;

    virtual void Stepper( const G4double yInput[],
                          const G4double dydx[],
                                G4double hstep,
                                G4double yOutput[],
                                G4double yError[] ) override;

    void Stepper( const G4double yInput[],
                  const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[],
                        G4double dydxOutput[] );

    virtual G4double DistChord() const override;
    virtual G4int IntegratorOrder() const override { return 4; }

  private:

    void makeStep( const G4double yInput[],
                   const G4double dydx[],
                   const G4double hstep,
                         G4double yOutput[],
                         G4double* dydxOutput = nullptr,
                         G4double* yError = nullptr) const;
  private:

    G4double fyIn[G4FieldTrack::ncompSVEC],
             fdydx[G4FieldTrack::ncompSVEC],
             fyOut[G4FieldTrack::ncompSVEC],
             fdydxOut[G4FieldTrack::ncompSVEC];

    G4double fhstep = -1.0;
};

#endif
