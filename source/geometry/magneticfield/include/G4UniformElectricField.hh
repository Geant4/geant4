// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UniformElectricField.hh,v 1.2 1999-12-15 14:49:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//   Class for creation of Uniform Magnetic Field
//
// 30.1.97 V.Grichine
//  1.8.97 J.Apostolakis, cleanup, new 3-vector constructor, 
//                        and removal of helix-stepper (to separate file)
// 5.11.97 G.Cosmo, added copy constructor and assignment operator.

#ifndef G4UNIFORMELECTRICFIELD_HH
#define G4UNIFORMELECTRICFIELD_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ElectroMagneticField.hh"
#include "G4Mag_EqRhs.hh"

class G4UniformElectricField : public G4ElectroMagneticField
{
public:
         //  A field with value equal to FieldVector
         // 
         G4UniformElectricField(const G4ThreeVector FieldVector );

         G4UniformElectricField(G4double vField,
	                   G4double vTheta,
		           G4double vPhi     ) ;
		       
         ~G4UniformElectricField() ;

         // Copy constructor and assignment operator
         //
         G4UniformElectricField(const G4UniformElectricField &p);
         G4UniformElectricField& operator = (const G4UniformElectricField &p);

         void GetFieldValue(const G4double position[] ,
	                          G4double B[]      ) const ;

protected:

private:
         G4double fFieldComponents[6] ;
} ;

#endif
