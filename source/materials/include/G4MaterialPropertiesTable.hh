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
//
////////////////////////////////////////////////////////////////////////
//
// class G4MaterialPropertiesTable
//
// Class description:
//
// A Material properties table is a hash table, with
// key = property name, and value either G4double or
// G4MaterialPropertyVector

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

#include "globals.hh"
#include "G4MaterialPropertiesIndex.hh"
#include "G4MaterialPropertyVector.hh"

#include <map>

class G4MaterialPropertiesTable
{
 public:
  G4MaterialPropertiesTable();
  virtual ~G4MaterialPropertiesTable();

  void AddConstProperty(const G4String& key, G4double PropertyValue,
                        G4bool createNewKey = false);
  void AddConstProperty(const char* key, G4double PropertyValue,
                        G4bool createNewKey = false);
  // Add a new property to the table by giving a key-name and value

  G4MaterialPropertyVector* AddProperty(
    const G4String& key, const std::vector<G4double>& photonEnergies,
    const std::vector<G4double>& propertyValues, G4bool createNewKey = false);
  // Add a new property to the table by giving a key-name and
  // vectors of values

  G4MaterialPropertyVector* AddProperty(const char* key,
                                        G4double* PhotonEnergies,
                                        G4double* PropertyValues,
                                        G4int NumEntries,
                                        G4bool createNewKey = false);
  // Add a new property to the table by giving a key-name and the
  // arrays x and y of size NumEntries.

  void AddProperty(const G4String& key, G4MaterialPropertyVector* opv,
                   G4bool createNewKey = false);
  void AddProperty(const char* key, G4MaterialPropertyVector* opv,
                   G4bool createNewKey = false);
  // Add a new property to the table by giving a key-name and an
  // already constructed G4MaterialPropertyVector.

  void AddProperty(const G4String& key, const G4String& mat);
  // Add a new property to the table by giving a key name and a material
  // name. Properties are in namespace G4OpticalMaterialProperties
  // Not possible to create a new key with this method.

  void RemoveConstProperty(const G4String& key);
  void RemoveConstProperty(const char* key);
  // Remove a constant property from the table.

  void RemoveProperty(const G4String& key);
  void RemoveProperty(const char* key);
  // Remove a property from the table.

  G4double GetConstProperty(const G4String& key) const;
  G4double GetConstProperty(const char* key) const;
  // Get the constant property from the table corresponding to the key-name

  G4double GetConstProperty(const G4int index) const;
  // Get the constant property from the table corresponding to the key-index

  G4bool ConstPropertyExists(const G4String& key) const;
  G4bool ConstPropertyExists(const char* key) const;
  // Return true if a const property 'key' exists.

  G4bool ConstPropertyExists(const G4int index) const;
  // Return true if a const property with key-index 'index' exists.

  G4MaterialPropertyVector* GetProperty(const char* key,
                                        G4bool warning = false);
  G4MaterialPropertyVector* GetProperty(const G4String& key,
                                        G4bool warning = false);
  // Get the property from the table corresponding to the key-name.

  G4MaterialPropertyVector* GetProperty(const G4int index,
                                        G4bool warning = false);
  // Get the property from the table corresponding to the key-index.

  void AddEntry(const G4String& key, G4double aPhotonEnergy,
                G4double aPropertyValue);
  void AddEntry(const char* key, G4double aPhotonEnergy,
                G4double aPropertyValue);
  // Add a new entry (pair of numbers) to the table for a given key.

  G4int GetConstPropertyIndex(const G4String& key,
                              G4bool warning = false) const;
  // Get the constant property index from the key-name

  G4int GetPropertyIndex(const G4String& key, G4bool warning = false) const;
  // Get the property index by the key-name.

  std::vector<G4String> GetMaterialPropertyNames() const;
  std::vector<G4String> GetMaterialConstPropertyNames() const;

  void DumpTable();

  const std::map<G4int, G4MaterialPropertyVector*, std::less<G4int>>*
  GetPropertyMap() const
  {
    return &MP;
  }
  const std::map<G4int, G4double, std::less<G4int>>* GetConstPropertyMap() const
  {
    return &MCP;
  }
  // Accessors required for persistency purposes

 private:
  G4MaterialPropertyVector* CalculateGROUPVEL();
  // Calculate the group velocity based on RINDEX

  std::map<G4int, G4MaterialPropertyVector*, std::less<G4int>> MP;
  typedef std::map<G4int, G4MaterialPropertyVector*,
                   std::less<G4int>>::const_iterator MPiterator;

  std::map<G4int, G4double, std::less<G4int>> MCP;
  typedef std::map<G4int, G4double, std::less<G4int>>::const_iterator
    MCPiterator;
  // material property map and constant property map by index types

  std::vector<G4String> G4MaterialPropertyName;
  std::vector<G4String> G4MaterialConstPropertyName;
  // vectors of strings of property names
};

#endif /* G4MaterialPropertiesTable_h */
