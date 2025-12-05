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
// G4UniformMagField
//
// Class description:
//
// Class for creation of Uniform Magnetic Field.

// Author: Vladimir Grichine (CERN), 30.01.1997
// -------------------------------------------------------------------
#ifndef G4UNIFORMMAGFIELD_HH
#define G4UNIFORMMAGFIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

/**
 * @brief G4UniformMagField is a class for defining a uniform magnetic field.
 */

class G4UniformMagField : public G4MagneticField
{
  public:
  
    /**
     * Constructor for G4UniformMagField, a field with value equal
     * to 'FieldVector'.
     *  @param[in] FieldVector The field vector value.
     */
    G4UniformMagField(const G4ThreeVector& FieldVector);

    /**
     * Alternative constructor for G4UniformMagField.
     *  @param[in] vField The field component.
     *  @param[in] vTheta The Theta component.
     *  @param[in] vPhi The Phi component.
     */
    G4UniformMagField(G4double vField,
                      G4double vTheta,
                      G4double vPhi);

    /**
     * Default Destructor.
     */
    ~G4UniformMagField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4UniformMagField(const G4UniformMagField& p);
    G4UniformMagField& operator = (const G4UniformMagField& p);

    /**
     * Returns the field value 'MagField' on given time 'yTrack'.
     */
    void GetFieldValue(const G4double yTrack[4],
                             G4double* MagField) const override;

    /**
     * Sets the field value.
     */
    void SetFieldValue(const G4ThreeVector& newFieldValue);

    /**
     * Returns the constant field value.
     */
    G4ThreeVector GetConstantFieldValue() const;
    
    /**
     * Returns a pointer to a new allocated clone of this object.
     */
    G4Field* Clone() const override;

  private:

    G4double fFieldComponents[3];
};

#endif
