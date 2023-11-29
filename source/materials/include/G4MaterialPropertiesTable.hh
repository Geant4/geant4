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

////////////////////////////////////////////////////////////////////////
//
// class G4MaterialPropertiesTable
//
// Class description:
//
// A Material properties table is a hash table, with
// key = property name, and value either G4double or
// G4MaterialPropertyVector
//
// File:        G4MaterialPropertiesTable.hh
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     2005-05-12 add SetGROUPVEL() by P. Gumplinger
//              2002-11-05 add named material constants by P. Gumplinger
//              1999-11-05 Migration from G4RWTPtrHashDictionary to STL
//                         by John Allison
//              1999-10-29 add method and class descriptors
//              1997-03-25 by Peter Gumplinger
//              > cosmetics (only)
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertiesTable_h
#define G4MaterialPropertiesTable_h 1

#include "G4MaterialPropertiesIndex.hh"
#include "G4MaterialPropertyVector.hh"
#include "globals.hh"

#include <vector>

class G4MaterialPropertiesTable
{
 public:
  G4MaterialPropertiesTable();
  virtual ~G4MaterialPropertiesTable();

  // Add a new property to the table by giving a key-name and value
  void AddConstProperty(const G4String& key, G4double propertyValue, G4bool createNewKey = false);
  void AddConstProperty(const char* key, G4double propertyValue, G4bool createNewKey = false);

  // Add a new property to the table by giving a key-name and
  // vectors of values
  G4MaterialPropertyVector* AddProperty(const G4String& key,
    const std::vector<G4double>& photonEnergies, const std::vector<G4double>& propertyValues,
    G4bool createNewKey = false, G4bool spline = false);

  // Add a new property to the table by giving a key-name and the
  // arrays x and y of size NumEntries.
  G4MaterialPropertyVector* AddProperty(const char* key, G4double* photonEnergies,
    G4double* propertyValues, G4int numEntries, G4bool createNewKey = false, G4bool spline = false);

  // Add a new property to the table by giving a key-name and an
  // already constructed G4MaterialPropertyVector.
  void AddProperty(const G4String& key, G4MaterialPropertyVector* opv, G4bool createNewKey = false);
  void AddProperty(const char* key, G4MaterialPropertyVector* opv, G4bool createNewKey = false);

  // Add a new property to the table by giving a key name and a material
  // name. Properties are in namespace G4OpticalMaterialProperties
  // Not possible to create a new key with this method.
  void AddProperty(const G4String& key, const G4String& mat);

  // Remove a constant property from the table.
  void RemoveConstProperty(const G4String& key);
  void RemoveConstProperty(const char* key);

  // Remove a property from the table.
  void RemoveProperty(const G4String& key);
  void RemoveProperty(const char* key);

  // Get a constant property from the table
  // It is an error to ask for a const property that the user has not defined.
  //  Check if it has been defined with ConstPropertyExists() first.
  G4double GetConstProperty(const G4String& key) const;
  G4double GetConstProperty(const char* key) const;
  G4double GetConstProperty(const G4int index) const;

  // Return true if a const property has been defined by the user.
  // Despite the name, this returns false for a const property in
  //  GetMaterialConstPropertyNames() but not defined by user.
  // Use this method before calling GetConstProperty().
  G4bool ConstPropertyExists(const G4String& key) const;
  G4bool ConstPropertyExists(const char* key) const;
  G4bool ConstPropertyExists(const G4int index) const;

  // Get the property from the table corresponding to the key-index or index.
  // nullptr is returned if the property has not been defined by the user.
  G4MaterialPropertyVector* GetProperty(const char* key) const;
  G4MaterialPropertyVector* GetProperty(const G4String& key) const;
  G4MaterialPropertyVector* GetProperty(const G4int index) const;

  // Add a new entry (pair of numbers) to the table for a given key.
  void AddEntry(const G4String& key, G4double aPhotonEnergy, G4double aPropertyValue);
  void AddEntry(const char* key, G4double aPhotonEnergy, G4double aPropertyValue);

  // Get the constant property index from the key-name
  // It is an error to request the index of a non-existent key (key not
  //  present in fMaterialConstPropertyNames()).
  G4int GetConstPropertyIndex(const G4String& key) const;

  // Get the property index by the key-name.
  // It is an error to request the index of a non-existent key (key not
  //  present in GetMaterialPropertyNames()).
  G4int GetPropertyIndex(const G4String& key) const;

  // print the material properties and material constant properties
  void DumpTable() const;

  // the next four methods are used in persistency/GDML:
  const std::vector<G4String>& GetMaterialPropertyNames() const { return fMatPropNames; }
  const std::vector<G4String>& GetMaterialConstPropertyNames() const { return fMatConstPropNames; }
  // return references to the vectors of material (constant) properties.
  const std::vector<G4MaterialPropertyVector*>& GetProperties() const { return fMP; }
  const std::vector<std::pair<G4double, G4bool>>& GetConstProperties() const { return fMCP; }

 private:
  // Calculate the group velocity based on RINDEX
  G4MaterialPropertyVector* CalculateGROUPVEL();

  // Vector of pointer to material property vectors.
  // All entries are initialized to nullptr.  Pointer is not null when mat.prop. vector defined.
  // Order of entries in MP defined by enum in G4MaterialPropertiesIndex.
  std::vector<G4MaterialPropertyVector*> fMP;

  // Vector of energy-independent (i.e., "constant") material properties. We
  // need to keep track if a property is defined or not: the bool in the pair
  // is 'true' if the property is defined.
  // Order of entries in MCP defined by enum in G4MaterialPropertiesIndex.
  std::vector<std::pair<G4double, G4bool>> fMCP;

  std::vector<G4String> fMatPropNames;  // vector of strings of property names
  std::vector<G4String> fMatConstPropNames;  // vector of strings of property names
};

#endif /* G4MaterialPropertiesTable_h */
