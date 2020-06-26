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
// G4VIsotopeTable
//
// Class description:
//
// Virtual class for access to an isotope table containing
// "stable" isotopes with their properties.

// Author: H.Kurashige, 5 October 1999
// --------------------------------------------------------------------
#ifndef G4VIsotopeTable_hh
#define G4VIsotopeTable_hh 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4Ions.hh"

class G4IsotopeProperty;

class G4VIsotopeTable
{
  public:

    G4VIsotopeTable();
    explicit G4VIsotopeTable(const G4String&);
      // Constructors

    G4VIsotopeTable(const G4VIsotopeTable&);
    G4VIsotopeTable& operator=(const G4VIsotopeTable&);
      // Copy contructor and assignment operator

    virtual ~G4VIsotopeTable();
      // Destructor

    virtual G4IsotopeProperty* GetIsotope(G4int Z, G4int A, G4double E,
            G4Ions::G4FloatLevelBase flb=G4Ions::G4FloatLevelBase::no_Float)= 0;
      // Pure virtual method

    virtual G4IsotopeProperty* GetIsotopeByIsoLvl(G4int Z, G4int A,
                                                  G4int level=0);
      // Search the isotope in the table. 
      // The isotope is designated by 
      //    G4int    Z:  number of proton (Atomic number)
      //    G4int    A:  number of nucleon (Atomic mass)
      //      and
      //    G4double E:  excited energy
      //    G4Ions::G4FloatLevelBase flb: floating level base (from G4Ions.hh)
      //      or
      //    G4int  level: isomer level
      // in the given G4IsotopeProperty.
      // If corresponding isotope exist in the table, this method returns
      // 'true', as well as fills other properties such as spin, lifetime,
      // decay modes and precise excited energy in the given G4IsotopeProperty.
      // This method returns 'false' if no corresponding isotope is found 
      // without modification of property

    G4int GetVerboseLevel() const;
    void SetVerboseLevel(G4int level);  
      // Set/Get verbose level

    void DumpTable(G4int Zmin=1, G4int Zmax=118);
      // Dump table

    const G4String& GetName() const;

  private:

    G4String fName = "";
    G4int verboseLevel = 0;
};

// ------------------------
// Inline methods
// ------------------------

inline
const G4String&  G4VIsotopeTable::GetName() const
{
  return fName;
}

inline 
G4int G4VIsotopeTable::GetVerboseLevel() const
{
  return verboseLevel;
}

inline 
void G4VIsotopeTable::SetVerboseLevel(G4int level)
{
  verboseLevel = level;
}

#endif
