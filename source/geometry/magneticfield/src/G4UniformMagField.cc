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
// $Id: G4UniformMagField.cc,v 1.10 2004/12/02 09:55:21 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// Class for creation of uniform Magnetic Field
//
// 30.1.97 V.Grichine
//
// -------------------------------------------------------------------

#include "G4UniformMagField.hh"

G4UniformMagField::G4UniformMagField(const G4ThreeVector& FieldVector )
{
      fFieldComponents[0] = FieldVector.x();
      fFieldComponents[1] = FieldVector.y();
      fFieldComponents[2] = FieldVector.z();
}

void
G4UniformMagField::SetFieldValue(const G4ThreeVector& newFieldVector )
{
      fFieldComponents[0] = newFieldVector.x();
      fFieldComponents[1] = newFieldVector.y();
      fFieldComponents[2] = newFieldVector.z();
}
   
G4UniformMagField::G4UniformMagField(G4double vField,
                                     G4double vTheta,
                                     G4double vPhi    )
{
   if(vField >= 0 && 
      vTheta >= 0 && vTheta <= pi && 
      vPhi >= 0 && vPhi <= twopi)
   {
      fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
      fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
      fFieldComponents[2] = vField*std::cos(vTheta) ;
   }
   else
   {
      G4Exception("G4UniformMagField::G4UniformMagField()",
                  "WrongArgumentValue", FatalException, "Invalid parameters.") ;
   }
}

G4UniformMagField::~G4UniformMagField()
{
}

G4UniformMagField::G4UniformMagField (const G4UniformMagField &p)
   : G4MagneticField(p)
{
   for (G4int i=0; i<3; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
}

G4UniformMagField& G4UniformMagField::operator = (const G4UniformMagField &p)
{
   if (&p == this) return *this;
   for (G4int i=0; i<3; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
   return *this;
}

// ------------------------------------------------------------------------

void G4UniformMagField::GetFieldValue (const G4double [4],
                                             G4double *B  ) const 
{
   B[0]= fFieldComponents[0] ;
   B[1]= fFieldComponents[1] ;
   B[2]= fFieldComponents[2] ;
}

G4ThreeVector G4UniformMagField::GetConstantFieldValue() const
{
   G4ThreeVector B(fFieldComponents[0],
                   fFieldComponents[1],
                   fFieldComponents[2]);
  return B;
}
