// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectroMagneticField.hh,v 1.2 1999-12-15 14:49:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A full Electromagnetic field, containing both electric and magnetic fields.
//
//   It is an abstract class, and a derived type of this field must be
//     created by the user to describe his/her field configuration.
//
//  Created: JA, November 12th, 1998
//                     

#ifndef G4ELECTROMAGNETIC_FIELD_DEF
#define G4ELECTROMAGNETIC_FIELD_DEF

#include "G4MagneticField.hh"

class G4ElectroMagneticField : public G4MagneticField
{
  public:

     G4ElectroMagneticField();
     virtual ~G4ElectroMagneticField();

     //  Copy constructor & assignment operator
     G4ElectroMagneticField(const G4ElectroMagneticField &p);
     G4ElectroMagneticField& operator = (const G4ElectroMagneticField &p);

     virtual void  GetFieldValue( const  double Point[3],
					 double *Bfield ) const = 0;
};

// Implementation 

inline  G4ElectroMagneticField::G4ElectroMagneticField() {}
inline  G4ElectroMagneticField::~G4ElectroMagneticField() {}
inline  G4ElectroMagneticField::G4ElectroMagneticField(const G4ElectroMagneticField &p) {}
 // Not needed: { *this = p; }

inline  G4ElectroMagneticField& 
G4ElectroMagneticField::operator = (const G4ElectroMagneticField &p)
 { *this = p; return *this; }

#endif /* G4ELECTROMAGNETIC_FIELD_DEF */
