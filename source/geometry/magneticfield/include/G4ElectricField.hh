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
// $Id: G4ElectricField.hh,v 1.1 2003/11/05 10:35:55 japost Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
// class G4ElectricField
//
// Class description:
//
// Electric Field abstract class, implements inquiry function interface.

// History:
// - Created. JA, November 4th, 2003.

#ifndef G4ELECTRIC_FIELD_DEF
#define G4ELECTRIC_FIELD_DEF

#include "G4Types.hh"
#include "G4ElectroMagneticField.hh"

class G4ElectricField : public G4ElectroMagneticField
{
  public:  // with description

     G4ElectricField();
     virtual ~G4ElectricField();
       // Constructor and destructor. No actions.

     G4ElectricField(const G4ElectricField &r);
     G4ElectricField& operator = (const G4ElectricField &p);
       // Copy constructor & assignment operator.

     G4bool   DoesFieldChangeEnergy() const { return true; }
       // Since an electric field can change track energy

     virtual void  GetFieldValue( const G4double Point[4],
					G4double *Bfield ) const = 0;
};

#endif /* G4ELECTRIC_FIELD_DEF */
