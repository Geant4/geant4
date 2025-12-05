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
// G4UniformElectricField
//
// Class description:
//
// Class for creation of Uniform electric Magnetic Field.

// Author: Vladimir Grichine (CERN), 30.01.1997
// -------------------------------------------------------------------
#ifndef G4UNIFORMELECTRICFIELD_HH
#define G4UNIFORMELECTRICFIELD_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4ElectricField.hh"

/**
 * @brief G4UniformElectricField is class defining a uniform
 * electric magnetic field.
 */

class G4UniformElectricField : public G4ElectricField
{
  public:

    /**
     * Constructor for G4UniformElectricField, a field with value equal
     * to 'FieldVector'.
     *  @param[in] FieldVector The field vector value.
     */
    G4UniformElectricField(const G4ThreeVector& FieldVector);

    /**
     * Alternative constructor for G4UniformElectricField.
     *  @param[in] vField The field component.
     *  @param[in] vTheta The Theta component.
     *  @param[in] vPhi The Phi component.
     */
    G4UniformElectricField(G4double vField,
                           G4double vTheta,
                           G4double vPhi);

    /**
     * Default Destructor.
     */
    ~G4UniformElectricField() override = default;

    /**
     * Copy constructor and assignment operator.
     */
    G4UniformElectricField(const G4UniformElectricField &p);
    G4UniformElectricField& operator = (const G4UniformElectricField &p);

    /**
     * Returns the field value 'field' on given time 'pos'.
     */
    void GetFieldValue(const G4double pos[4], G4double* field) const override;

    /**
     * Returns a pointer to a new allocated clone of this object.
     */
    G4Field* Clone() const override;

  private:
  
    G4double fFieldComponents[6];
};

#endif
