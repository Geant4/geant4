// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02ElectroMagneticField.hh,v 1.2 2001-04-04 15:11:15 grichine Exp $
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


#ifndef F02ElectroMagneticField_H
#define F02ElectroMagneticField_H

#include "G4MagneticField.hh"
// #include "G4ElectroMagneticField.hh"
#include "G4UniformElectricField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class F02FieldMessenger;

class F02ElectroMagneticField :  public G4UniformElectricField
{
public:
  F02ElectroMagneticField(G4ThreeVector) ;  //  The value of the field
  F02ElectroMagneticField() ;               //  A zero field

 ~F02ElectroMagneticField() ;  
      
  void  GetFieldValue( const  G4double Point[3],
			      G4double Bfield[3] ) const {;} ;
  
  void SetStepperType( G4int i) { fStepperType = i ; } ;

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; } ;

  void UpdateField();

  void SetFieldValue(G4ThreeVector fieldVector) ;
  void SetFieldValue(G4double      fieldValue) ;
  G4ThreeVector GetConstantFieldValue();

protected:

      // Find the global Field Manager

  G4FieldManager*         GetGlobalFieldManager() ;   // static 

  G4FieldManager*         fFieldManager ;

  G4ChordFinder*          fChordFinder ;

  G4EqMagElectricField*   fEquation ;
 
  G4UniformElectricField* fEMfield ; 

  G4MagIntegratorStepper* fStepper ;

  G4int                   fStepperType ;

  G4double                fMinStep ;
 
  F02FieldMessenger*      fFieldMessenger;

};

#endif
