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
// Class for creation of Uniform Gravitation Field.
//

// History:
// - 14.06.11 P.Gumplinger, Created.
// -------------------------------------------------------------------
// Adopted from G4UniformElectricField.hh
//
// Thanks to Peter Fierlinger (PSI) and
// A. Capra and A. Fontana (INFN Pavia)
// -------------------------------------------------------------------
//
#include "G4UniformGravityField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Construct from a 3-vector
G4UniformGravityField::G4UniformGravityField(const G4ThreeVector FieldVector)
  : G4Field ( true ) //  Gravity flag *on*
{
      fFieldComponents[0] = FieldVector.x();
      fFieldComponents[1] = FieldVector.y();
      fFieldComponents[2] = FieldVector.z();
}

// Construct from a double > default = -9.81 m*s^-2
G4UniformGravityField::G4UniformGravityField(const G4double gy )
  : G4Field ( true ) 
{
      fFieldComponents[0] = 0.0;
      fFieldComponents[1] = gy;
      fFieldComponents[2] = 0.0;
}

G4Field* G4UniformGravityField::Clone() const
{
    return new G4UniformGravityField( G4ThreeVector(fFieldComponents[0],
                                                    fFieldComponents[1],
                                                    fFieldComponents[2]) );
}

G4UniformGravityField::~G4UniformGravityField()
{
}

G4UniformGravityField::G4UniformGravityField (const G4UniformGravityField &p)
   : G4Field(p)
{
   for (G4int i=0; i<3; i++)
   {
     fFieldComponents[i] = p.fFieldComponents[i];
   }
}

G4UniformGravityField&
G4UniformGravityField::operator = (const G4UniformGravityField &p)
{
  if (&p == this) return *this; 
  G4Field::operator=(p); 
  for (G4int i=0; i<3; i++)
  {
    fFieldComponents[i] = p.fFieldComponents[i];
  }
  return *this;
}

// -------------------------------------------------------------------

void G4UniformGravityField::GetFieldValue (const G4double [4],
                                                 G4double *G ) const
{
   G[0]= fFieldComponents[0] ;
   G[1]= fFieldComponents[1] ;
   G[2]= fFieldComponents[2] ;
}
