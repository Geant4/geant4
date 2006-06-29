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
// $Id: F03ElectroMagneticField.hh,v 1.2 2006-06-29 18:29:16 gunter Exp $
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


#ifndef F03ElectroMagneticField_H
#define F03ElectroMagneticField_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class F03FieldMessenger;

class F03ElectroMagneticField: public G4MagneticField
{
public:
  F03ElectroMagneticField(G4ThreeVector) ;  //  The value of the field
  F03ElectroMagneticField() ;               //  A zero field

 ~F03ElectroMagneticField() ;  
      
  void  GetFieldValue( const  G4double Point[3],
			      G4double *Bfield ) const {;} ;
  
  void SetStepperType( G4int i) { fStepperType = i ; } ;

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; } ;

  void UpdateField();

  void SetFieldValue(G4ThreeVector fieldVector) ;
  void SetFieldValue(G4double      fieldValue) ;
  G4ThreeVector GetConstantFieldValue();
  G4FieldManager*  GetLocalFieldManager() { return fLocalFieldManager ;};
protected:

      // Find the global Field Manager

  G4FieldManager*         GetGlobalFieldManager() ;   // static 

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
