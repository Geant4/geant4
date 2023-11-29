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
// G4Mag_EqRhs
//
// Class description:
//
// The "standard" equation of motion of a particle in a pure magnetic field.
// Other that might be required are:
//     i) is when using a moving reference frame ... or
//    ii) extending for other forces, e.g. an electric field

// Created: J.Apostolakis, CERN - 13.01.1997
// --------------------------------------------------------------------
#ifndef G4MAG_EQRHS_HH
#define G4MAG_EQRHS_HH

#include "G4Types.hh"
#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

class G4MagneticField;

class G4Mag_EqRhs : public G4EquationOfMotion
{
  public:

    G4Mag_EqRhs(G4MagneticField* magField);
   ~G4Mag_EqRhs() override;
      // Constructor and destructor. No actions.

    void EvaluateRhsGivenB( const G4double y[],
                            const G4double B[3],
                                  G4double dydx[] ) const override = 0;
      // Given the value of the  field "B", this function 
      // calculates the value of the derivative dydx.
      // This is the _only_ function a subclass must define.
      // The other two functions use Rhs_givenB.

    inline G4double FCof() const { return fCof_val; }

    void SetChargeMomentumMass( G4ChargeState particleCharge,
                                G4double MomentumXc,
                                G4double mass ) override;
  private:

    G4double fCof_val = 0.0;

    static const G4double fUnitConstant;     // Set to 0.299792458
      // Coefficient in the Lorentz motion equation (Lorentz force), if the
      // magnetic field B is in Tesla, the particle charge in units of the 
      // elementary (positron?) charge, the momentum P in MeV/c, and the
      // space coordinates and path along the trajectory in mm.
};

#endif
