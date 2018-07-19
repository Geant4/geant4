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
// $Id: G4UniformElectricField.cc 96751 2016-05-04 09:39:38Z gcosmo $
//
// 
//
// Class for creation of uniform Electric Field
//
// 30.1.97 V.Grichine
//
// -------------------------------------------------------------------

#include "G4UniformElectricField.hh"
#include "G4PhysicalConstants.hh"

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
   if ( (vField<0) || (vTheta<0) || (vTheta>pi) || (vPhi<0) || (vPhi>twopi) )
   {
      G4Exception("G4UniformElectricField::G4UniformElectricField()",
                  "GeomField0002", FatalException, "Invalid parameters.");
   }

   fFieldComponents[0] = 0.0;
   fFieldComponents[1] = 0.0;
   fFieldComponents[2] = 0.0;
   fFieldComponents[3] = vField*std::sin(vTheta)*std::cos(vPhi) ;
   fFieldComponents[4] = vField*std::sin(vTheta)*std::sin(vPhi) ;
   fFieldComponents[5] = vField*std::cos(vTheta) ;
}

G4Field* G4UniformElectricField::Clone() const
{
    return new G4UniformElectricField( G4ThreeVector(fFieldComponents[3],
                                                     fFieldComponents[4],
                                                     fFieldComponents[5]) );
}

G4UniformElectricField::~G4UniformElectricField()
{
}

G4UniformElectricField::G4UniformElectricField (const G4UniformElectricField &p)
   : G4ElectricField(p)
{
   for (G4int i=0; i<6; i++)
   {
     fFieldComponents[i] = p.fFieldComponents[i];
   }
}

G4UniformElectricField&
G4UniformElectricField::operator = (const G4UniformElectricField &p)
{
  if (&p == this) return *this; 
  G4ElectricField::operator=(p); 
  for (G4int i=0; i<6; i++)
  {
    fFieldComponents[i] = p.fFieldComponents[i];
  }
  return *this;
}

// ------------------------------------------------------------------------

void G4UniformElectricField::GetFieldValue (const G4double[4],
                                            G4double *fieldBandE ) const 
{
   fieldBandE[0]= 0.0;
   fieldBandE[1]= 0.0;
   fieldBandE[2]= 0.0;
   fieldBandE[3]= fFieldComponents[3] ;
   fieldBandE[4]= fFieldComponents[4] ;
   fieldBandE[5]= fFieldComponents[5] ;
}
