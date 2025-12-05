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
// G4UniformGravifyField
//
// Class description:
//
// Class for creation of Uniform Gravitation Field.

// Author: Peter Gumplinger (TRIUMF), 14.06.2011
//         Thanks to P.Fierlinger (PSI), A.Capra and A.Fontana (INFN Pavia)
// -------------------------------------------------------------------
#ifndef G4UNIFORMGRAVITYFIELD_HH
#define G4UNIFORMGRAVITYFIELD_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Field.hh"

/**
 * @brief G4UniformGravityField is a class for defining a uniform
 * gravitation field.
 */

class G4UniformGravityField : public G4Field
{
  public:

    /**
     * Constructor for G4UniformGravityField, a field with value equal
     * to 'FieldVector'.
     *  @param[in] FieldVector The field vector value.
     */
    G4UniformGravityField(const G4ThreeVector& FieldVector );

    /**
     * Alternative constructor for G4UniformGravityField.
     *  @param[in] gy The gravitation field value (default is the Standard
     *             Gravitational field on earth's surface.
     */
    G4UniformGravityField(const G4double gy = -9.81*CLHEP::m/CLHEP::s/CLHEP::s);

    /**
     * Default Destructor.
     */
    ~G4UniformGravityField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4UniformGravityField(const G4UniformGravityField& p);
    G4UniformGravityField& operator=(const G4UniformGravityField& p);

    /**
     * The field can change track energy, so returning true.
     */
    inline G4bool DoesFieldChangeEnergy() const override { return true; }

    /**
     * Returns the field value 'field' on given time 'Point'.
     */
    void GetFieldValue(const G4double Point[4], G4double* field) const override;
    
    /**
     * Returns a pointer to a new allocated clone of this object.
     */
    G4Field* Clone() const override;

    /**
     * Returns the field type ID, 'kGravity'.
     */
    inline G4FieldType GetFieldType() const override { return kGravity; }

  private:

    G4double fFieldComponents[3];
};

#endif
