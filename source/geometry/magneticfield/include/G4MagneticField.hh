//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4MagneticField.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
