// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Mag_UsualEqRhs.hh,v 1.1 1999-01-07 16:07:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  This is the standard right-hand side for equation of motion.
//
//    The only case another is required is when using a moving reference
//     frame ... or extending the class to include additional Forces,
//     eg an electric field
//
//            J. Apostolakis, January 13th, 1997
//
#ifndef G4MAG_USUAL_EQRHS
#define G4MAG_USUAL_EQRHS

#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"

class G4Mag_UsualEqRhs: public G4Mag_EqRhs{
   public:
     G4Mag_UsualEqRhs( G4MagneticField* MagField ) :
		     G4Mag_EqRhs( MagField ) {};
    ~G4Mag_UsualEqRhs() {} ;  

     //  Given the value of the magnetic field B, this function 
     //   calculates the value of the derivative dydx.
     //
     void EvaluateRhsGivenB( const  G4double y[],
			      const  G4double B[3],
				     G4double dydx[] ) const;
};

#endif /* G4MAG_USUAL_EQRHS */
