// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagneticField.hh,v 1.2 1999-12-15 14:49:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Magnetic Field abstract class,  implements inquiry function interface.
//
//       JA, January 13th, 1996
//
//  November 5th, 1997 - G.Cosmo, added default & copy constructors, virtual
//                       destructor and assignment operator.

#ifndef G4MAGNETIC_FIELD_DEF
#define G4MAGNETIC_FIELD_DEF

#include "G4Field.hh"

class G4MagneticField : public G4Field
{
  public:

     G4MagneticField();
     virtual ~G4MagneticField();

     //  Copy constructor & assignment operator
     G4MagneticField(const G4MagneticField &p);
     G4MagneticField& operator = (const G4MagneticField &p);

     //  Old version of field evaluation function:
     //  to be replaced by following function (GetFieldValue)
     // virtual void MagneticField( const  double Point[3],
     //					double Bfield[3] ) = 0;

     virtual void  GetFieldValue( const  double Point[3],
					 double *Bfield ) const = 0;
};

// Implementation 

inline  G4MagneticField::G4MagneticField() {}
inline  G4MagneticField::~G4MagneticField() {}
inline  G4MagneticField::G4MagneticField(const G4MagneticField &p) {}
 // Not needed: { *this = p; }
inline  G4MagneticField& G4MagneticField::operator = (const G4MagneticField &p)
 { *this = p; return *this; }

#endif /* G4MAGNETIC_FIELD_DEF */
