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

//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2015 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file G4FieldBuilder.h
/// \brief Definition of the G4FieldBuilder class
///
/// \author I. Hrivnacova; IJCLab, Orsay

#ifndef G4FIELDBUILDER_HH
#define G4FIELDBUILDER_HH

#include "G4Cache.hh"
#include "G4FieldParameters.hh"
#include "globals.hh"

#include <vector>

class G4FieldBuilderMessenger;
class G4FieldSetup;
class G4LogicalVolume;
class G4EquationOfMotion;
class G4MagIntegratorStepper;

/// \brief The manger class for building magnetic or other field
/// using the configuration in field parameters.
///
/// Purpose: Provide a single 'place' to configure field & integration
///
/// - It can configure a global field, and field(s) local to a (logical) volume
/// - The parameter values can be configured by the user (else use a default)
/// - They can be set/changed via a messenger provided or in the code
///      of the user detector construciton
/// - It retains ownership of the following object(s):
///      field parameters and field setups, field
///
///  Note MT: an object of the builder class should be created on master only
///      (in DetectorConstruction constructor)
///      The functions SetGlobal/LocalField and ConstructFieldSetup should
///      be called on workers (in DetectorConstruction::ConstructSDandField )
///
///  This design/implementation covers the most common use cases.
///  It cannot be used to create some complex setups such as
///    - equations templated on the field type,
///    - steppers/drivers templated on the equation and field types.
///
/// \author I. Hrivnacova; IJCLab, Orsay

class G4FieldBuilder
{
 public:
  /// Destructor
  ~G4FieldBuilder();

  // Static access method
  //

  /// Create the class instance, if it does not exist,
  /// and return it on the next calls.
  static G4FieldBuilder* Instance();

  /// Return the information if an instance exists
  static G4bool IsInstance();

  // Functions for constructing field setup
  //

  /// Create local magnetic field parameters (configuration) which can be then
  /// configured by the user via UI commands.
  /// The parameters are used in geometry only if a local magnetic field is
  /// associated with the volumes with the given name
  G4FieldParameters* CreateFieldParameters(const G4String& fieldVolName);

  /// Construct setups for all registered fields.
  void ConstructFieldSetup();

  /// Update magnetic field.
  /// This function must be called if the field parameters were changed
  /// in other than PreInit> phase.
  void UpdateField();

  /// Reinitialize if geometry has been modified.
  /// This function is called by G4RunManager during ReinitializeGeometry()
  void Reinitialize();

  // Set methods
  //

  /// Default field type is set to kMagnetic;
  /// this function should be called for other than magnetic field
  /// in order to update the default equation and stepper types.
  void SetFieldType(G4FieldType fieldType);

  // Set or reset the global field.
  // Update field objects, if the field was already constructed.
  // If warn, issue a warning if the previous field is deleted.
  void SetGlobalField(G4Field* field, G4bool warn = false);

  /// Register the local field in the map.
  /// Update field objects, if the field was already constructed.
  /// If warn, issue a warning if the previous field is deleted.
  /// The field is propagated to all volume daughters regardless
  /// if they have already assigned a field manager or not.
  /// When multiple local fields are defined (by calling this function
  /// multiple times), they will be applied in the order they were set.
  void SetLocalField(G4Field* field, G4LogicalVolume* lv, G4bool warn = false);

  /// Set user equation of motion
  void SetUserEquationOfMotion(
    G4EquationOfMotion* equation, G4String volumeName = "");

  /// Set user stepper
  void SetUserStepper(
    G4MagIntegratorStepper* stepper, G4String volumeName = "");

  /// Set verbose level
  void SetVerboseLevel(G4int value);

  // Get methods
  //

  /// Get field parameters with the given volumeName.
  /// Return global field parameters, if volume name is empty.
  G4FieldParameters* GetFieldParameters(const G4String& volumeName = "") const;

 private:
  /// Default constructor
  G4FieldBuilder();
  /// Not implemented
  G4FieldBuilder(const G4FieldBuilder& right) = delete;
  /// Not implemented
  G4FieldBuilder& operator=(const G4FieldBuilder& right) = delete;

  // Methods

  /// Get field parameters with the given volumeName or create them if they
  /// do not exist yet
  G4FieldParameters* GetOrCreateFieldParameters(const G4String& volumeName);

  /// Get field setup with the given logical volume
  G4FieldSetup* GetFieldSetup(G4LogicalVolume* lv);

  /// Create magnetic, electromagnetic or gravity field setup
  void CreateFieldSetup(G4Field* field,
    G4FieldParameters* fieldParameters, G4LogicalVolume* lv);

  /// Construct Geant4 global magnetic field setup
  void ConstructGlobalField();

  /// Construct Geant4 local magnetic field setups from the local fields map
  void ConstructLocalFields();

  /// Update all field setups
  void UpdateFieldSetups();

  // helper methods
  std::vector<G4FieldSetup*>& GetFieldSetups();
  std::vector<std::pair<G4LogicalVolume*, G4Field*>>& GetLocalFields();

  // Data members

  /// Information if an instance exists
  inline static G4ThreadLocal G4bool fgIsInstance { false };

  /// Messenger for this class
  G4FieldBuilderMessenger* fMessenger = nullptr;

  /// Field parameters
  std::vector<G4FieldParameters*> fFieldParameters;

  /// Field setups
  G4Cache<std::vector<G4FieldSetup*>*> fFieldSetups;

  /// Registered global field
  static G4ThreadLocal G4Field* fGlobalField;

  /// Registered local fields
  G4Cache<std::vector<std::pair<G4LogicalVolume*, G4Field*>>*> fLocalFields;

  /// info if field objects were constructed
  static G4ThreadLocal G4bool fIsConstructed;

  /// verbose level
  G4int fVerboseLevel = 1;
};

// inline methods

inline G4bool G4FieldBuilder::IsInstance()
{
  // Return the information if an instance exists
  return fgIsInstance;
}

inline void G4FieldBuilder::SetVerboseLevel(G4int value)
{
  // Set verbose level
  fVerboseLevel = value;
}

inline std::vector<G4FieldSetup*>& G4FieldBuilder::GetFieldSetups()
{
  // Return reference to field setups from G4Cache
  return *fFieldSetups.Get();
}

inline std::vector<std::pair<G4LogicalVolume*, G4Field*>>& G4FieldBuilder::GetLocalFields()
{
  // Return reference to local fields map from G4Cache
  return *fLocalFields.Get();
}

#endif // G4FIELDBUILDER_HH
