// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Mag_EqRhs.hh,v 1.2 1999-02-12 12:29:42 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// The right hand size of the equation of motion of a particle 
//   in a magnetic field.
//
//   (Possible use of alternative to "normal" version: rotating reference 
//    frame)
//
//   JA, January 13th, 1996
//
#ifndef G4_MAG_EQRHS_DEF
#define G4_MAG_EQRHS_DEF
#include  "globals.hh"

#include  "G4EquationOfMotion.hh"

#include "G4MagneticField.hh"      // class G4MagneticField;  not enough ??

class G4Mag_EqRhs : public G4EquationOfMotion
{
  public:
     G4Mag_EqRhs( G4MagneticField *magField );
    ~G4Mag_EqRhs();

     //  Given the value of the  field "B", this function 
     //   calculates the value of the derivative dydx.
     //  --------------------------------------------------------
     //  This is the _only_ function a subclass must define.
     //  The other two functions use Rhs_givenB.
     //
     virtual void EvaluateRhsGivenB( const  G4double y[],
			      const  G4double B[3],
				     G4double dydx[] ) const = 0;

     G4double FCof() const { return fCof_val; }

     virtual void  SetChargeMomentumMass( const G4double particleCharge, // in e+ units
			          const G4double MomentumXc,
				  const G4double mass);
     
  private:

     G4double        fCof_val;

     // Coefficient in the Lorentz motion equation (Lorentz force), if the
     //  magnetic field B is in Tesla, the particle charge in units of the 
     //  elementary (positron?) charge, the momentum P in MeV/c, and the
     //  space coordinates and path along the trajectory in mm .
     //
     static const G4double fUnitConstant;     // Set in G4Mag_EqRhs.cc 
					      // to 0.299792458
};

#endif /* G4_MAG_EQRHS_DEF */
