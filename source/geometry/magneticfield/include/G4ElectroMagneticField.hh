// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectroMagneticField.hh,v 1.4 2000-11-01 15:15:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4ElectroMagneticField
//
// Class description:
//
// A full Electromagnetic field, containing both electric and magnetic fields.
// It is an abstract class, and a derived type of this field must be
// created by the user to describe his/her field configuration.

// History:
// - Created: J.Apostolakis, November 12th, 1998                   

#ifndef G4ELECTROMAGNETIC_FIELD_DEF
#define G4ELECTROMAGNETIC_FIELD_DEF

#include "G4MagneticField.hh"

class G4ElectroMagneticField : public G4MagneticField
{
  public:  // with description

     G4ElectroMagneticField() {;}
     virtual ~G4ElectroMagneticField() {;}

     G4ElectroMagneticField(const G4ElectroMagneticField &) {;}
     G4ElectroMagneticField& operator = (const G4ElectroMagneticField &);
       // Copy constructor & assignment operators.

     virtual void  GetFieldValue(const G4double Point[3],
				       G4double *Bfield ) const = 0;
};

// Inline implementations

inline
G4ElectroMagneticField& 
G4ElectroMagneticField::operator = (const G4ElectroMagneticField &p)
{
  if (&p == this) return *this;
  *this = p; return *this;
}

#endif /* G4ELECTROMAGNETIC_FIELD_DEF */
