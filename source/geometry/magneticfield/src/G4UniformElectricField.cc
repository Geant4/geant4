// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UniformElectricField.cc,v 1.3 2001-03-27 09:43:22 japost Exp $
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
      fFieldComponents[0] = 0.0;
      fFieldComponents[1] = 0.0;
      fFieldComponents[2] = 0.0;
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
      fFieldComponents[0] = 0.0;
      fFieldComponents[1] = 0.0;
      fFieldComponents[2] = 0.0;
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
   for (G4int i=0; i<6; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
}

G4UniformElectricField& G4UniformElectricField::operator = (const G4UniformElectricField &p)
{
   for (G4int i=0; i<6; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
   return *this;
}

// ------------------------------------------------------------------------


void G4UniformElectricField::GetFieldValue (const G4double [3],
                                             G4double fieldBandE[6]  ) const 
{
   fieldBandE[0]= 0.0;
   fieldBandE[1]= 0.0;
   fieldBandE[2]= 0.0;
   fieldBandE[3]= fFieldComponents[3] ;
   fieldBandE[4]= fFieldComponents[4] ;
   fieldBandE[5]= fFieldComponents[5] ;
   return ;
}

