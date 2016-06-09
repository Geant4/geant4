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
// $Id: G4MagneticField.hh,v 1.13 2003/11/05 10:35:55 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// class G4MagneticField
//
// Class description:
//
// Magnetic Field abstract class, implements inquiry function interface.

// History:
// - Created. JA, January 13th, 1996.
// --------------------------------------------------------------------

#ifndef G4MAGNETIC_FIELD_DEF
#define G4MAGNETIC_FIELD_DEF

#include "G4Types.hh"
#include "G4ElectroMagneticField.hh"

class G4MagneticField : public G4ElectroMagneticField
{
  public:  // with description

     G4MagneticField();
     virtual ~G4MagneticField();
       // Constructor and destructor. No actions.

     G4MagneticField(const G4MagneticField &r);
     G4MagneticField& operator = (const G4MagneticField &p);
       // Copy constructor & assignment operator.

     G4bool   DoesFieldChangeEnergy() const { return false; }
       //  Since a pure magnetic field does not change track energy

     virtual void  GetFieldValue( const G4double Point[4],
                                        G4double *Bfield ) const = 0;
};

#endif /* G4MAGNETIC_FIELD_DEF */
