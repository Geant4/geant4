//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ElectroMagneticField.hh,v 1.10 2004/09/01 10:00:10 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// class G4ElectroMagneticField
//
// Class description:
//
// A full Electromagnetic field, containing both electric and magnetic fields.
// It is an abstract class, and a derived type of this field must be
// created by the user to describe his/her field configuration.

// We have established a convention for the electromagnetic field components:
// In the GetValue method, the return values of Bfield will have 
// the following meaning
//  - Components 0, 1 and 2 are the Magnetic Field (x, y, z respectively);
//  - Components 3, 4 and 5 are the Electric field (x, y, z respectively).
// 
// Note 1: one or the other field could optional, depending on the Equation
// Note 2: such a convention is required between any field and its 
// corresponding equation of motion.
//
// History:
// - Created:  J. Apostolakis, November 12th, 1998
// - Modified: 
//   V. Grichine     8 Nov 2001: Extended "Point" arg to [4] array to add time
//   G. Cosmo        2 Apr 2003: Reorgansised, moved inline methods to .cc
//   J. Apostolakis  5 Nov 2003: Derive directly from G4Field 
//   J. Apostolakis 31 Aug 2004: Information on convention for components
// -------------------------------------------------------------------

#ifndef G4ELECTROMAGNETIC_FIELD_DEF
#define G4ELECTROMAGNETIC_FIELD_DEF

#include "G4Field.hh"

class G4ElectroMagneticField : public G4Field
{
  public:  // with description

     G4ElectroMagneticField();
     virtual ~G4ElectroMagneticField();

     G4ElectroMagneticField(const G4ElectroMagneticField &r);
     G4ElectroMagneticField& operator = (const G4ElectroMagneticField &p);
       // Copy constructor & assignment operators.

     virtual void  GetFieldValue(const G4double Point[4],
                                       G4double *Bfield ) const = 0;
       //  Return as Bfield[0], [1], [2] the magnetic field x, y & z components
       //   and   as Bfield[3], [4], [5] the electric field x, y & z components

     virtual G4bool   DoesFieldChangeEnergy() const = 0;

       //  For field with an electric component this should be true
       //  For pure magnetic field this should be false
       //    Alternative: default safe implementation { return true; }
};

#endif /* G4ELECTROMAGNETIC_FIELD_DEF */
