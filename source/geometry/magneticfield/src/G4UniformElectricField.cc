// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UniformElectricField.cc,v 1.2 1999-12-15 14:49:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Class for creation of uniform Electric Field
//
// 30.1.97 V.Grichine
//
#include "G4UniformElectricField.hh"
#include "globals.hh"
#include "geomdefs.hh"

G4UniformElectricField::G4UniformElectricField(const G4ThreeVector FieldVector )
{
      fFieldComponents[3] = FieldVector.x();
      fFieldComponents[4] = FieldVector.y();
      fFieldComponents[5] = FieldVector.z();
}
   
G4UniformElectricField::G4UniformElectricField(G4double vField,
			             G4double vTheta,
			             G4double vPhi    )
{
   if(vField >= 0 && 
      vTheta >= 0 && vTheta <= pi && 
      vPhi >= 0 && vPhi <= twopi)
   {
      fFieldComponents[3] = vField*sin(vTheta)*cos(vPhi) ;
      fFieldComponents[4] = vField*sin(vTheta)*sin(vPhi) ;
      fFieldComponents[5] = vField*cos(vTheta) ;
   }
   else
   {
      G4Exception("Invalid parameters in G4UniformElectricField::G4UniformElectricField") ;
   }
}

G4UniformElectricField::~G4UniformElectricField()
{
   ;
}

G4UniformElectricField::G4UniformElectricField (const G4UniformElectricField &p)
{
   for (G4int i=3; i<6; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
}

G4UniformElectricField& G4UniformElectricField::operator = (const G4UniformElectricField &p)
{
   for (G4int i=3; i<6; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
   return *this;
}

// ------------------------------------------------------------------------


void G4UniformElectricField::GetFieldValue (const G4double [3],
                                             G4double E[3]  ) const 
{
   E[0]= fFieldComponents[3] ;
   E[1]= fFieldComponents[4] ;
   E[2]= fFieldComponents[5] ;
   return ;
}

