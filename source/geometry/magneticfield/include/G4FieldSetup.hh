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

/// \file G4FieldSetup.h
/// \brief Definition of the G4FieldSetup class
///
/// This code was initially developed in Geant4 VMC package
/// (https://github.com/vmc-project)
/// and adapted to Geant4.
///
/// \author I. Hrivnacova; IJCLab, Orsay

#ifndef G4FIELDSETUP_HH
#define G4FIELDSETUP_HH

#include "G4FieldParameters.hh"
#include "globals.hh"

class G4Field;
class G4FieldParameters;
class G4FieldSetupMessenger;

class G4ChordFinder;
class G4EquationOfMotion;
class G4FieldManager;
class G4MagIntegratorStepper;
class G4LogicalVolume;
class G4VIntegrationDriver;

class TVirtualMagField;

/// \ingroup geometry
/// \brief The class for constructing magnetic, electromagnetic and gravity
/// fields which strength is defined via G4Field.
///
/// The equation of motion of a particle in a field and the
/// integration method is set according to the selection in
/// G4FieldParameters, as well as other accuracy parameters.
/// The default values in G4FieldParameters correspond to defaults
/// set in Geant4 (taken from Geant4 9.3 release.)
/// As Geant4 classes to not provide access methods for these defaults,
/// the defaults have to be checked with each new Geant4 release.
/// TO DO: unify defaults in G4 classes and G4 parameters
///
/// \author I. Hrivnacova; IJClab, Orsay

class G4FieldSetup
{
 public:
  /// Standard constructor
  G4FieldSetup(const G4FieldParameters& parameters, G4Field* field,
    G4LogicalVolume* lv = nullptr);
  /// Destructor
  ~G4FieldSetup();

  // Methods

  /// Clear previously created setup
  void Clear();
  /// Update field setup with new field parameters
  void Update();
  /// Print information
  void PrintInfo(G4int verboseLevel, const G4String about = "created");

  // Set methods

  /// Set G4 field
  void SetG4Field(G4Field* field);

  // Access to field setting

  /// Return the instantiated field
  G4Field* GetG4Field() const;
  /// Return the logical vol;ume
  G4LogicalVolume* GetLogicalVolume() const;
  /// Return the equation of motion
  G4EquationOfMotion* GetEquation() const;
  /// Return the magnetic integrator stepper
  G4MagIntegratorStepper* GetStepper() const;
  /// Return the magnetic integrator driver
  G4VIntegrationDriver* GetIntegrationDriver() const;

 private:
  /// Not implemented
  G4FieldSetup() = delete;
  /// Not implemented
  G4FieldSetup(const G4FieldSetup& right) = delete;
  /// Not implemented
  G4FieldSetup& operator=(const G4FieldSetup& right) = delete;

  // Methods

  // Create cached magnetic field if const distance is set > 0.
  // and field is of G4MagneticField.
  // Return the input field otherwise.
  G4Field* CreateCachedField(
    const G4FieldParameters& parameters, G4Field* field);

  /// Set the equation of motion of a particle in a field
  G4EquationOfMotion* CreateEquation(G4EquationType equation);

  /// Set the integrator of particle's equation of motion
  G4MagIntegratorStepper* CreateStepper(
    G4EquationOfMotion* equation, G4StepperType stepper);

  /// Set the FSAL integrator of particle's equation of motion
  G4VIntegrationDriver* CreateFSALStepperAndDriver(
    G4EquationOfMotion* equation, G4StepperType stepper, G4double minStep);

  // methods to update field setup step by step
  /// Create cached field (if ConstDistance is set)
  void CreateCachedField();
  /// Create cached field (if ConstDistance is set)
  void CreateStepper();
  /// Create chord finder
  void CreateChordFinder();
  /// Update field manager
  void UpdateFieldManager();

  // Data members

  /// Messenger for this class
  G4FieldSetupMessenger* fMessenger = nullptr;
  /// Parameters
  const G4FieldParameters& fParameters;
  /// Geant4 field manager
  G4FieldManager* fFieldManager = nullptr;
  /// Geant4 field
  G4Field* fG4Field = nullptr;
  /// The associated ROOT volume (if local field)
  G4LogicalVolume* fLogicalVolume = nullptr;
  /// The equation of motion
  G4EquationOfMotion* fEquation = nullptr;
  /// The magnetic integrator stepper
  G4MagIntegratorStepper* fStepper = nullptr;
  /// The magnetic integrator driver
  G4VIntegrationDriver* fDriver = nullptr;
  /// Chord finder
  G4ChordFinder* fChordFinder = nullptr;
};

// inline functions

inline void G4FieldSetup::SetG4Field(G4Field* field)
{
  // Set G4 field
  fG4Field = field;
}

inline G4Field* G4FieldSetup::GetG4Field() const
{
  // Return the instantiated field
  return fG4Field;
}

inline G4LogicalVolume* G4FieldSetup::GetLogicalVolume() const
{
  // Return the logical vol;ume
  return fLogicalVolume;
}

inline G4EquationOfMotion* G4FieldSetup::GetEquation() const
{
  // Return the equation of motion
  return fEquation;
}

inline G4MagIntegratorStepper* G4FieldSetup::GetStepper() const
{
  // Return the magnetic integrator stepper
  return fStepper;
}

#endif // G4FIELDSETUP_HH
