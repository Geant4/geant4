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
// G4ElectricField
//
// Class description:
//
// Electric field abstract class, implements inquiry function interface.

// Author: John Apostolakis (CERN), 04.11.2003
// --------------------------------------------------------------------
#ifndef G4ELECTRIC_FIELD_HH
#define G4ELECTRIC_FIELD_HH

#include "G4Types.hh"
#include "G4ElectroMagneticField.hh"

/**
 * @brief G4ElectricField is an abstract class for electric field.
 * It implements inquiry function interface.
 */

class G4ElectricField : public G4ElectroMagneticField
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4ElectricField() = default;
    ~G4ElectricField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4ElectricField(const G4ElectricField& r) = default;
    G4ElectricField& operator = (const G4ElectricField& p);

    /**
     * Returns true, since an electric field can change track energy.
     */
    inline G4bool DoesFieldChangeEnergy() const override { return true; }

    /**
     * Interface for returning the field value 'Bfield' on given time 'Point'.
     */
    void GetFieldValue( const G4double Point[4],
                              G4double* Bfield ) const override = 0;
};

#endif
