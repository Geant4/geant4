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
// $Id: G4ElectroMagneticField.hh,v 1.9 2003/11/05 10:37:06 japost Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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

     virtual G4bool   DoesFieldChangeEnergy() const = 0;

       //  For field with an electric component this should be true
       //  For pure magnetic field this should be false
       //    Alternative: default safe implementation { return true; }
};

#endif /* G4ELECTROMAGNETIC_FIELD_DEF */
