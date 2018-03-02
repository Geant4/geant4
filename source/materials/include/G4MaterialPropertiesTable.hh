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
// $Id: G4MaterialPropertiesTable.hh 108517 2018-02-16 08:18:55Z gcosmo $
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
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertiesTable_h
#define G4MaterialPropertiesTable_h 1

/////////////
// Includes
/////////////

#include <cmath>
#include <map>
#include "globals.hh"
#include "G4MaterialPropertiesIndex.hh"
#include "G4MaterialPropertyVector.hh"

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertiesTable
{
  public: // Without description

    G4MaterialPropertiesTable(); 
    virtual ~G4MaterialPropertiesTable();

  public: // With description

    inline void AddConstProperty(const char     *key,
                                 G4double PropertyValue);
    // Add a new property to the table by giving a key-name and value 

    G4MaterialPropertyVector* AddProperty(const char     *key,
                                          G4double *PhotonEnergies,
                                          G4double *PropertyValues,
                                          G4int     NumEntries);
    // Add a new property to the table by giving a key-name and the
    // arrays x and y of size NumEntries.

    void AddProperty(const char *key, G4MaterialPropertyVector *opv);
    // Add a new property to the table by giving a key-name and an
    // already constructed G4MaterialPropertyVector.

    inline void RemoveConstProperty(const char *key);
    // Remove a constant property from the table.

    inline void RemoveProperty(const char *key);
    // Remove a property from the table.

    G4double GetConstProperty(const char *key) const;
    // Get the constant property from the table corresponding to the key-name

    G4double GetConstProperty(const G4int index) const;
    // Get the constant property from the table corresponding to the key-index

    G4bool ConstPropertyExists(const char *key) const;
    // Return true if a const property 'key' exists.

    G4MaterialPropertyVector* GetProperty(const char *key,
                                          G4bool warning=false);
    // Get the property from the table corresponding to the key-name.

    G4MaterialPropertyVector* GetProperty(const G4int index, 
                                          G4bool warning=false);
    // Get the property from the table corresponding to the key-index.

    void AddEntry(const char *key, G4double aPhotonEnergy,
                                   G4double aPropertyValue);
    // Add a new entry (pair of numbers) to the table for a given key.

    G4int GetConstPropertyIndex(const G4String& key,
                                G4bool warning=false) const;
    // Get the constant property index from the key-name

    G4int GetPropertyIndex(const G4String& key, G4bool warning=false) const;
    // Get the property index by the key-name.

    std::vector<G4String> GetMaterialPropertyNames() const;
    std::vector<G4String> GetMaterialConstPropertyNames() const;

    void DumpTable();

  public:  // without description

    std::map< G4String, G4MaterialPropertyVector*, std::less<G4String> >*
      GetPropertiesMap(); 
    std::map< G4String, G4double, std::less<G4String> >*
      GetPropertiesCMap();

    const std::map<G4int, G4MaterialPropertyVector*, std::less<G4int> >*
      GetPropertyMap() const { return &MP; }
    const std::map<G4int, G4double, std::less<G4int> >*
      GetConstPropertyMap() const { return &MCP; }
    // Accessors required for persistency purposes

  private:

    G4MaterialPropertyVector* CalculateGROUPVEL();
    // Calculate the group velocity based on RINDEX

    G4MaterialPropertyVector* SetGROUPVEL();
    // Dummy method: will be obsolete from the next (version 11) release

  private:

    std::map<G4String, G4MaterialPropertyVector*, std::less<G4String> > MPT;
    typedef std::map< G4String, G4MaterialPropertyVector*,
                      std::less<G4String> >::const_iterator MPTiterator;

    std::map< G4String, G4double, std::less<G4String> > MPTC;
    typedef std::map< G4String, G4double,
                      std::less<G4String> >::const_iterator MPTCiterator;
    // MPT and MPTC will be obsolete when associate public interfaces, 
    // GetPropertiesMap and GetPropertiesCMap are removed from the version 11. 

    std::map<G4int, G4MaterialPropertyVector*, std::less<G4int> > MP;
    typedef std::map< G4int, G4MaterialPropertyVector*,
                      std::less<G4int> >::const_iterator MPiterator;

    std::map< G4int, G4double, std::less<G4int> > MCP;
    typedef std::map< G4int, G4double,
                      std::less<G4int> >::const_iterator MCPiterator;
    //material property map and constant property map by index types

    std::vector<G4String> G4MaterialPropertyName;
    std::vector<G4String> G4MaterialConstPropertyName;
    // vectors of strings of property names

};

/////////////////////
// Inline definitions
/////////////////////

#include "G4MaterialPropertiesTable.icc"

#endif /* G4MaterialPropertiesTable_h */
