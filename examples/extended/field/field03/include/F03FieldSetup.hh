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
/// \file field/field03/include/F03FieldSetup.hh
/// \brief Definition of the F03FieldSetup class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F03FieldSetup_h
#define F03FieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class F03FieldMessenger;

///  A class for setting up the Magnetic Field
///
///  It also creates the necessary classes to control accuracy of propagation.
///  In this example
///    - There is a global field for most of the setup;
///    - A local field overides it for some volume(s) and it assumed to be
///      uniform.

class F03FieldSetup
{
public:
  F03FieldSetup();           //  A zero field
  virtual ~F03FieldSetup();

  void SetStepperType( G4int i ) { fStepperType = i; }

  void CreateSteppers();

  void SetMinStep(G4double s) { fMinStep = s; }

  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldZValue(G4double      fieldValue);
  void SetLocalFieldValue(G4ThreeVector fieldVector);
  G4ThreeVector GetGlobalFieldValue() const { return GetConstantFieldValue(fMagneticField); }
  G4ThreeVector GetLocalFieldValue() const { return GetConstantFieldValue(fLocalMagneticField); }

  void UpdateField();

  G4FieldManager* GetLocalFieldManager() { return fLocalFieldManager;}

protected:

  // Find the global Field Manager
  G4FieldManager*         GetGlobalFieldManager() ;
  G4ThreeVector           GetConstantFieldValue(G4MagneticField* magneticField) const;

  G4FieldManager*         fFieldManager;
  G4FieldManager*         fLocalFieldManager;
  G4ChordFinder*          fChordFinder;
  G4ChordFinder*          fLocalChordFinder;
  G4Mag_UsualEqRhs*       fEquation;
  G4Mag_UsualEqRhs*       fLocalEquation;
  G4MagneticField*        fMagneticField;
  G4MagneticField*        fLocalMagneticField;

  G4MagIntegratorStepper* fStepper;
  G4MagIntegratorStepper* fLocalStepper;
  G4int                   fStepperType;

  G4double                fMinStep;
 
  F03FieldMessenger*      fFieldMessenger;

};

#endif
