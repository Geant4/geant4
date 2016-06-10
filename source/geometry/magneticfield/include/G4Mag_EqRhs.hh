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
// $Id: G4Mag_EqRhs.hh 69699 2013-05-13 08:50:30Z gcosmo $
//
//
// class G4Mag_EqRhs
//
// Class description:
//
// The "standard" equation of motion of a particle in a pure magnetic field.

// History:
// - Created. J.Apostolakis, January 13th 1997
// --------------------------------------------------------------------

#ifndef G4_MAG_EQRHS_DEF
#define G4_MAG_EQRHS_DEF

#include "G4Types.hh"
#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

class G4MagneticField;

class G4Mag_EqRhs : public G4EquationOfMotion
{
  public: // with description

     G4Mag_EqRhs( G4MagneticField *magField );
     virtual ~G4Mag_EqRhs();
       // Constructor and destructor. No actions.

     virtual void EvaluateRhsGivenB( const  G4double y[],
                                     const  G4double B[3],
                                            G4double dydx[] ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

     inline G4double FCof() const;

     virtual void SetChargeMomentumMass( G4ChargeState particleCharge,
                                         G4double MomentumXc,
                                         G4double mass);
     
  private:

     G4double fCof_val;

     static const G4double fUnitConstant;     // Set in G4Mag_EqRhs.cc 
                                              // to 0.299792458
       // Coefficient in the Lorentz motion equation (Lorentz force), if the
       //  magnetic field B is in Tesla, the particle charge in units of the 
       //  elementary (positron?) charge, the momentum P in MeV/c, and the
       //  space coordinates and path along the trajectory in mm .
};

inline
G4double G4Mag_EqRhs::FCof() const
{
  return fCof_val;
}

#endif /* G4_MAG_EQRHS_DEF */
