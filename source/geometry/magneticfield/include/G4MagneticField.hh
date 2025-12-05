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
// G4MagneticField
//
// Class description:
//
// Magnetic Field abstract class, implements inquiry function interface.

// Author: John Apostolakis (CERN), 13.01.1996
// --------------------------------------------------------------------
#ifndef G4MAGNETIC_FIELD_HH
#define G4MAGNETIC_FIELD_HH

#include "G4Types.hh"
#include "G4Field.hh"

class G4MagneticField : public G4Field
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4MagneticField();
    ~G4MagneticField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4MagneticField(const G4MagneticField& r);
    G4MagneticField& operator= (const G4MagneticField& p);

    /**
     * Since a pure magnetic field does not change track energy, returns false.
     */
    inline G4bool DoesFieldChangeEnergy() const override { return false; }

    /**
     * Given the position time vector 'Point', returns the value of the
     * field in the array 'Bfield'.
     *  @param[in] Point The position time vector.
     *  @param[out] Bfield The field array in output.
     */
    void GetFieldValue( const G4double Point[4],
                              G4double* Bfield ) const override = 0;

    /**
     * Returns the field type-ID, "kMagnetic".
     * This should be overriden in derived classes.
     */
    inline G4FieldType GetFieldType() const override { return kMagnetic; }
};

#endif
