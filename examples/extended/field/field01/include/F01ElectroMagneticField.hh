// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01ElectroMagneticField.hh,v 1.1 2001-03-27 16:21:26 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
// 
//  $ Id:  $

// Should this be a:
//    i) messenger
//   ii) user class that creates the field       ? 
//  iii) simply a derived class of Uniform field ?  <== I have chosen this now.
//   iv) a field manager that creates/updates field    (Prefered?)
//


#ifndef F01ElectroMagneticField_H
#define F01ElectroMagneticField_H

#include "G4UniformMagField.hh"

class G4FieldManager;

class F01ElectroMagneticField: public G4UniformMagField
{
public:
      F01ElectroMagneticField(G4ThreeVector) ;  //  The value of the field
      F01ElectroMagneticField() ;               //  A zero field

     ~F01ElectroMagneticField() ;  
      
      // Set the field to (0, 0, fieldValue)

      void SetFieldValue(G4ThreeVector fieldVector) ;
      void SetFieldValue(G4double      fieldValue) ;
      G4ThreeVector GetConstantFieldValue();

protected:

      // Find the global Field Manager

      G4FieldManager* GetGlobalFieldManager() ;   // static 
      G4int fStepperType ;
};

#endif
