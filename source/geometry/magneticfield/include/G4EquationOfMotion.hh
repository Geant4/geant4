// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EquationOfMotion.hh,v 1.2 1999-12-15 14:49:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Abstract Base Class for:
//
//   The right hand size of the equation of motion of a particle in a field.
//
#ifndef G4_EquationOfMotion_DEF
#define G4_EquationOfMotion_DEF
#include  "globals.hh"

#include  "G4Field.hh"    // class G4Field;  // should be enough

class G4EquationOfMotion 
{
  public:
     G4EquationOfMotion( G4Field *Field );
     virtual ~G4EquationOfMotion();

     //  Given the value of the  field "B", this function 
     //   calculates the value of the derivative dydx.
     //  --------------------------------------------------------
     //  This is the _only_ function a subclass must define.
     //  The other two functions use Rhs_givenB.
     //
     virtual void EvaluateRhsGivenB( const  G4double y[],
			      const  G4double B[3],
				     G4double dydx[] ) const = 0;

     //  Set the charge, momentum and mass of the current particle
     //   --> used to set the equation's coefficients ...
     virtual void  SetChargeMomentumMass( 
				 const G4double particleCharge, // in e+ units
			         const G4double MomentumXc,
			         const G4double MassXc2) = 0;
     

     //  This calculates the value of the derivative dydx at y.
     //  It is the usual enquiry function.
     //  ---------------------------
     //  (It is not virtual, but calls the virtual function above.)
     //
     void RightHandSide( const  G4double y[],
				G4double dydx[]   ) const;

     //  Same as RHS above, but also returns the value of B.
     //
     //   Should be made the new default ? after putting dydx & B in a class
     //
     void EvaluateRhsReturnB( const  G4double y[],
			      G4double dydx[],
			      G4double Field[]  ) const;

     //  Obtain only the field - the stepper assumes it is pure Magnetic
     //     Not protected, because G4RKG3_Stepper uses it directly
     
     void GetFieldValue( const  G4double Point[3],
			        G4double Field[] )  const
			   {  itsField-> GetFieldValue( Point, Field ); }

     G4Field* GetFieldObj();
     void     SetFieldObj(G4Field* pField);

  //------------------------------------------------------------------------
  //public:
  //  virtual void doNothing(); // To help compiler with their virtual tables.

  private:

     G4Field *itsField;

};

#include "G4EquationOfMotion.icc"

#endif /* G4_EquationOfMotion_DEF */
