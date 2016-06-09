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
// $Id: F02ElectricFieldSetup.hh,v 1.2 2006-06-29 17:17:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
