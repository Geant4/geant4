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
//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2015 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file G4FieldBuilder.cc
/// \brief Implementation of the G4FieldBuilder class
///
/// \author I. Hrivnacova; IJCLab, Orsay

#include "G4FieldBuilder.hh"
#include "G4FieldBuilderMessenger.hh"

#include "G4Exception.hh"
#include "G4FieldParameters.hh"
#include "G4FieldSetup.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4ThreadLocal G4Field*  G4FieldBuilder::fGlobalField = nullptr;
G4ThreadLocal G4bool G4FieldBuilder::fIsConstructed = false;

//_____________________________________________________________________________
G4FieldBuilder* G4FieldBuilder::Instance()
{
  static G4FieldBuilder instance;
  fgIsInstance = true;
  return &instance;
}

//_____________________________________________________________________________
G4FieldBuilder::G4FieldBuilder()
{
  // Default constructor
 
  fMessenger = new G4FieldBuilderMessenger(this);

  // Field parameters for global field
  fFieldParameters.push_back(new G4FieldParameters());
}

//_____________________________________________________________________________
G4FieldBuilder::~G4FieldBuilder()
{
  // Destructor

  delete fMessenger;

  for (auto parameters : fFieldParameters) {
    delete parameters;
  }

  for (auto setup : GetFieldSetups()) {
    delete setup;
  }

  fgIsInstance = false;

  // magnetic field objects are deleted via G4 kernel
}

//
// private methods
//

//_____________________________________________________________________________
G4FieldParameters* G4FieldBuilder::GetOrCreateFieldParameters(
  const G4String& volumeName)
{
  // Get field parameters with the given volumeName or create them if they
  // do not exist yet

  // Get user field parameters
  for (auto fieldParameters : fFieldParameters) {
    if (fieldParameters->GetVolumeName() == volumeName) {
      // G4cout << "Found field parameters for " << volumeName << G4endl;
      return fieldParameters;
    }
  }

  // Create field parameters if not yet defined
  // G4cout << "Go to create field parameters for " << volumeName << G4endl;
  auto fieldParameters = new G4FieldParameters(volumeName);
  fFieldParameters.push_back(fieldParameters);
  return fieldParameters;
}

//_____________________________________________________________________________
G4FieldSetup* G4FieldBuilder::GetFieldSetup(G4LogicalVolume* lv)
{
  // Get field setup with the given logical volume

  // Get user field parameters
  for (auto fieldSetup : GetFieldSetups()) {
    if (fieldSetup->GetLogicalVolume() == lv) {
      return fieldSetup;
    }
  }

  return nullptr;
}

//_____________________________________________________________________________
void G4FieldBuilder::CreateFieldSetup(G4Field* field,
  G4FieldParameters* fieldParameters, G4LogicalVolume* lv)
{
  // Create magnetic, electromagnetic or gravity field setup

  auto fieldSetup = new G4FieldSetup(*fieldParameters, field, lv);
  fieldSetup->PrintInfo(fVerboseLevel);

  if (fFieldSetups.Get() == nullptr) {
    auto fieldSetups = new std::vector<G4FieldSetup*>();
    fFieldSetups.Put(fieldSetups);
  }

  GetFieldSetups().push_back(fieldSetup);
}

//_____________________________________________________________________________
void G4FieldBuilder::ConstructGlobalField()
{
  // Construct Geant4 global magnetic field setup

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::ConstructGlobalField " << G4endl;
  }

  CreateFieldSetup(fGlobalField, fFieldParameters[0], nullptr);
}

//_____________________________________________________________________________
void G4FieldBuilder::ConstructLocalFields()
{
  // Construct Geant4 local magnetic field setups from the local fields map

  if (fLocalFields.Get() == nullptr) return;

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::ConstructLocalFields()" << G4endl;
  }

  // Loop over local field map
  for (auto [lv, field] : GetLocalFields()) {

    // Volume name
    const auto& volumeName = lv->GetName();

    // Get or create user field parameters
    G4FieldParameters* fieldParameters =
      GetOrCreateFieldParameters(volumeName);

    if (fVerboseLevel > 1) {
      G4cout << "Construct local field in volume: " << volumeName << G4endl;
    }

    // Create magnetic field
    CreateFieldSetup(field, fieldParameters, lv);
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::UpdateFieldSetups()
{
  // Update all field setups

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::UpdateFieldSetups " << G4endl;
  }

  for (auto fieldSetup : GetFieldSetups()) {
    fieldSetup->Update();

    if (fVerboseLevel > 1) {
      fieldSetup->PrintInfo(fVerboseLevel, "updated");
    }
  }
}

//
// public methods
//

//_____________________________________________________________________________
G4FieldParameters* G4FieldBuilder::CreateFieldParameters(
  const G4String& fieldVolName)
{
  // Create local magnetic field parameters (configuration) which can be then
  // configured by the user via UI commands.
  // The parameters are used in geometry only if a local magnetic field is
  // associated with the volumes with the given name.

  auto fieldParameters = new G4FieldParameters(fieldVolName);
  fFieldParameters.push_back(fieldParameters);

  return fieldParameters;
}

//_____________________________________________________________________________
void G4FieldBuilder::ConstructFieldSetup()
{
  // Construct setups for all registered fields.

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::ConstructField" << G4endl;
  }

  if (fIsConstructed) {
    G4Exception(
      "G4FieldBuilder::ConstructField:", "GeomFieldParameters0001",
      JustWarning, "Field was already constructed.");
    return;
  }

  ConstructGlobalField();
  ConstructLocalFields();

  UpdateFieldSetups();

  fIsConstructed = true;
}

//_____________________________________________________________________________
void G4FieldBuilder::UpdateField()
{
  // Update magnetic field.
  // This function must be called if the field parameters were changed
  // in other than PreInit> phase.

  if (fFieldSetups.Get() == nullptr) {
    G4Exception(
      "G4FieldBuilder::UpdateField", "GeomFieldParameters0001",
      JustWarning, "No field setup is defined.");
    return;
  }

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::UpdateField" << G4endl;
  }

  // Update the objects defined in G4FieldSetup's
  UpdateFieldSetups();
}

//_____________________________________________________________________________
void G4FieldBuilder::Reinitialize()
{
  // Reinitialize if geometry has been modified.
  // This function is called by G4RunManager during ReinitializeGeometry()

  if (fVerboseLevel > 1) {
    G4cout << "G4FieldBuilder::Reinitialize" << G4endl;
  }
  
  // Delete global field
  delete fGlobalField;
  fGlobalField = nullptr;

  // Delete local fields if defined
  if (fLocalFields.Get() != nullptr) {
    for (auto vectorElement : GetLocalFields()) {
      delete vectorElement.second;
    }
    // Clear local fields map
    GetLocalFields().clear();
  }

  // Clear field setups if defined
  if (fFieldSetups.Get() != nullptr) {
    for (auto fieldSetup : GetFieldSetups()) {
      fieldSetup->SetG4Field(nullptr);
      fieldSetup->Clear();
    }
  }

  fIsConstructed = false;

  if (fVerboseLevel > 1) {
    G4cout << "End of G4FieldBuilder::Reinitialize" << G4endl;
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::SetFieldType(G4FieldType fieldType)
{
// Default field type is set to kMagnetic;
// this function should be called for other than magnetic field
// in order to update the default equation and stepper types.

  if (fIsConstructed) {
    G4Exception(
      "G4FieldBuilder::SetFieldType:", "GeomFieldParameters0001",
      JustWarning, "Field was already constructed.");
    return;
  }

  fFieldParameters[0]->SetFieldType(fieldType);

  // change default equation and stepper if other than magnetic field
  if (fieldType == kElectroMagnetic) {
    fFieldParameters[0]->SetEquationType(kEqElectroMagnetic);
    fFieldParameters[0]->SetStepperType(kClassicalRK4);
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::SetGlobalField(G4Field* field, G4bool warn)
{
  // Set or reset the global field.
  // Update field objects, if the field was already constructed.
  // If warn, issue a warning if the previous field is deleted.

  if (fGlobalField != nullptr && warn) {
    G4Exception(
      "G4FieldBuilder::SetGlobalField:", "GeomFieldParameters0001",
      JustWarning, "The global field already exists, it will be deleted.");
  }
  delete fGlobalField;
  fGlobalField = field;

  if (fIsConstructed) {
    // update the global field objects if already constructed
    GetFieldSetups()[0]->SetG4Field(field);
    GetFieldSetups()[0]->Update();
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::SetLocalField(
  G4Field* field, G4LogicalVolume* lv, G4bool warn)
{
  // Register the local field in the map.
  // Update field objects, if the field was already constructed.
  // If warn, issue a warning if the previous field is deleted.

  if (lv == nullptr) {
    G4cerr << "Cannot register local field without Logical volume." << G4endl;
    return;
  }

  if (fLocalFields.Get() == nullptr) {
    auto localFields = new std::vector<std::pair<G4LogicalVolume*, G4Field*>>();
    fLocalFields.Put(localFields);
  }

  auto it = GetLocalFields().begin();
  for (it = GetLocalFields().begin(); it != GetLocalFields().end(); ++it) {
    if (it->first == lv) break;
  }

  if (it != GetLocalFields().end()) {
    // replaced field if already in the map
    if (warn) {
      G4ExceptionDescription descr;
      descr << "Logical volume " << lv->GetName() << " has already field."
        " It will be deleted.";
      G4Exception(
        "G4FieldBuilder::SetLocalField:", "GeomFieldParameters0001",
        JustWarning, descr);
    }
    delete it->second;
    it->second = field;
  }
  else {
    // register field in the map
    GetLocalFields().push_back(std::pair(lv,field));
  }

  if (fIsConstructed) {
    // update this local field objects if already constructed
    auto fieldSetup = GetFieldSetup(lv);
    if (fieldSetup == nullptr) {
      G4Exception(
        "G4FieldBuilder::SetLocalField:", "GeomFieldParameters0001",
        JustWarning, "Cannot get field setup for a local field.");
      return;
    }
    fieldSetup->SetG4Field(field);
    fieldSetup->Update();
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::SetUserEquationOfMotion(
  G4EquationOfMotion* equation, G4String volumeName)
{
  // Set user equation of motion

  if (!volumeName.size()) {
    // global field
    fFieldParameters[0]->SetUserEquationOfMotion(equation);
  }
  else {
    // local field
    // Get or create user field parameters
    G4FieldParameters* fieldParameters =
      GetOrCreateFieldParameters(volumeName);
    fieldParameters->SetUserEquationOfMotion(equation);
  }
}

//_____________________________________________________________________________
void G4FieldBuilder::SetUserStepper(
  G4MagIntegratorStepper* stepper, G4String volumeName)
{
  // Set user stepper

  if (!volumeName.size()) {
    // global field
    fFieldParameters[0]->SetUserStepper(stepper);
  }
  else {
    // local field
    // Get or create user field parameters
    G4FieldParameters* fieldParameters =
      GetOrCreateFieldParameters(volumeName);
    fieldParameters->SetUserStepper(stepper);
  }
}

//_____________________________________________________________________________
G4FieldParameters* G4FieldBuilder::GetFieldParameters(
  const G4String& volumeName) const
{
  // Get field parameters with the given volumeName.
  // Return global field parameters, if volume name is empty.

  // Get user field parameters
  for (auto fieldParameters : fFieldParameters) {
    if (fieldParameters->GetVolumeName() == volumeName) {
      // G4cout << "Found field parameters for " << volumeName << G4endl;
      return fieldParameters;
    }
  }

  G4Exception(
    "G4FieldBuilder::GetFieldParameters:", "GeomFieldParameters0001",
    JustWarning, "Field parameters not found.");
  return nullptr;
}
