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
// $Id: F03ElectroMagneticField.hh,v 1.3 2001-10-15 17:20:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
// 

#ifndef F03ElectroMagneticField_H
#define F03ElectroMagneticField_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class F03FieldMessenger;

class F03ElectroMagneticField : public G4MagneticField
{

public:

  F03ElectroMagneticField() ;               //  A zero field
  F03ElectroMagneticField(G4ThreeVector) ;  //  The value of the field
 ~F03ElectroMagneticField() ;
      
  void  GetFieldValue( const  G4double Point[3],
			      G4double *Bfield ) const {;}
  
  void SetStepperType( G4int i) { fStepperType = i ; }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; }

  void UpdateField();

  void SetFieldValue(G4ThreeVector fieldVector) ;
  void SetFieldValue(G4double      fieldValue) ;
  G4ThreeVector GetConstantFieldValue();
  G4FieldManager*  GetLocalFieldManager() { return fLocalFieldManager ;}

protected:

  G4FieldManager*         GetGlobalFieldManager() ;
    // Returns the global Field Manager

  G4FieldManager*         fFieldManager ;
  G4FieldManager*         fLocalFieldManager ;

  G4ChordFinder*          fChordFinder ;
  G4ChordFinder*          fLocalChordFinder ;

  G4Mag_UsualEqRhs*       fEquation ; 
  G4Mag_UsualEqRhs*       fLocalEquation ; 

  G4MagneticField*        fMagneticField ; 
  G4MagneticField*        fLocalMagneticField ; 

  G4MagIntegratorStepper* fStepper ;
  G4MagIntegratorStepper* fLocalStepper ;
  G4int                   fStepperType ;

  G4double                fMinStep ;
 
  F03FieldMessenger*      fFieldMessenger;

};

#endif
