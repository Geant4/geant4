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
// $Id: F02ElectricFieldSetup.hh,v 1.1 2003/11/12 16:36:14 japost Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//    A class for control of the Electric Field of the detector.
//  The field for this case is uniform.
//
//  It is simply a 'setup' class that creates the field and necessary other parts
//


#ifndef F02ElectricFieldSetup_H
#define F02ElectricFieldSetup_H

#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver; 
class F02FieldMessenger;
class F02ElectricFieldSetup 
{
public:

  F02ElectricFieldSetup(G4ThreeVector) ;  //  The value of the field
  F02ElectricFieldSetup() ;               //  A zero field

 ~F02ElectricFieldSetup() ;  
      
  void SetStepperType( G4int i) { fStepperType = i ; }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; }

  void UpdateField();

  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();

protected:

      // Find the global Field Manager
  G4FieldManager*         GetGlobalFieldManager() ;

private:
  G4FieldManager*         fFieldManager ;

  G4ChordFinder*          fChordFinder ;

  G4EqMagElectricField*   fEquation ;

  G4ElectricField*        fEMfield;
 
  G4ThreeVector           fElFieldValue ; 

  G4MagIntegratorStepper* fStepper ;
  G4MagInt_Driver*        fIntgrDriver;

  G4int                   fStepperType ;

  G4double                fMinStep ;
 
  F02FieldMessenger*      fFieldMessenger;

};

#endif
