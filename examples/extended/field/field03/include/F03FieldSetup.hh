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
// $Id: F03FieldSetup.hh,v 1.1 2003/12/01 17:03:02 japost Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//  A class for setting up the Magnetic Field of the setup, and 
//   creating the necessary classes to control accuracy of propagation.
//  In this example
//    - There is a global field for most of the setup;
//    - A local  field overides it for some volume(s) and it assumed to be uniform.
// 

#ifndef F03FieldSetup_H
#define F03FieldSetup_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class F03FieldMessenger;

class F03FieldSetup
{

public:

  F03FieldSetup() ;               //  A zero field
  F03FieldSetup(G4ThreeVector) ;  //  The value of the field
 ~F03FieldSetup() ;
      
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
