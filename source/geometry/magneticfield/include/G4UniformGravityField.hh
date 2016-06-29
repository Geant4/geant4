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
// class G4UniformGravifyField
//
// Class description:
//
// Class for creation of Uniform Gravitation Field.
//

// History:
// - 14.06.11 P.Gumplinger, Created.
// -------------------------------------------------------------------
// Adapted from G4UniformElectricField.hh
//
// Thanks to Peter Fierlinger (PSI) and
// A. Capra and A. Fontana (INFN Pavia)
// -------------------------------------------------------------------
//
#ifndef G4UNIFORMGRAVITYFIELD_HH
#define G4UNIFORMGRAVITYFIELD_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Field.hh"

class G4UniformGravityField : public G4Field
{
  public:  // with description

    G4UniformGravityField(const G4ThreeVector FieldVector );
      // A field with value equal to FieldVector.

    G4UniformGravityField(const G4double gy = -9.81*CLHEP::m/CLHEP::s/CLHEP::s);
      // Standard Gravitational field on earth's surface

    virtual ~G4UniformGravityField();

    G4UniformGravityField(const G4UniformGravityField &p);
    G4UniformGravityField& operator = (const G4UniformGravityField &p);
      // Copy constructor and assignment operator

    G4bool   DoesFieldChangeEnergy() const { return true; }
      // Since a gravitational field can change track energy

    virtual void GetFieldValue(const G4double Point[4], G4double *field) const;
    
    virtual G4Field* Clone() const;

  private:

    G4double fFieldComponents[3];
};

#endif
