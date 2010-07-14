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
// $Id: G4UniformMagField.cc,v 1.12 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
   if ( (vField<0) || (vTheta<0) || (vTheta>pi) || (vPhi<0) || (vPhi>twopi) )
   {
      G4Exception("G4UniformMagField::G4UniformMagField()",
                  "WrongArgumentValue", FatalException, "Invalid parameters.") ;
   }
   fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
   fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
   fFieldComponents[2] = vField*std::cos(vTheta) ;
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
