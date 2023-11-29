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
// G4UnitsTable
//
// Class description:
//
// This class maintains a table of Units.
// A Unit has a name, a symbol, a value and belong to a category (i.e. its
// dimensional definition): Length, Time, Energy, etc...
// The Units are grouped by category. The TableOfUnits is a list of categories.
// The class G4BestUnit allows to convert automaticaly a physical quantity
// from its internal value into the most appropriate Unit of the same category.

// Author: M.Maire, 17.05.1998 - First version
// Revisions: G.Cosmo, 06.03.2001 - Migrated to STL vectors
// --------------------------------------------------------------------
#ifndef G4UnitsTable_hh
#define G4UnitsTable_hh 1

#include <vector>

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4UnitsCategory;
class G4UnitDefinition;

// --------------------------------------------------------------------

#ifdef G4MULTITHREADED

class G4UnitsTable : public std::vector<G4UnitsCategory*>
{
 public:
  using std::vector<G4UnitsCategory*>::vector;
  G4UnitsTable() = default;
  ~G4UnitsTable();

 public:
  void Synchronize();
  G4bool Contains(const G4UnitDefinition*, const G4String&);
};
#else

using G4UnitsTable = std::vector<G4UnitsCategory*>;

#endif

// --------------------------------------------------------------------

class G4UnitDefinition
{
 public:
  G4UnitDefinition(const G4String& name, const G4String& symbol,
                   const G4String& category, G4double value);

  ~G4UnitDefinition() = default;
  G4bool operator==(const G4UnitDefinition&) const;
  G4bool operator!=(const G4UnitDefinition&) const;

  inline const G4String& GetName() const;
  inline const G4String& GetSymbol() const;
  inline G4double GetValue() const;

  void PrintDefinition();

  static void BuildUnitsTable();
  static void PrintUnitsTable();
  static void ClearUnitsTable();

  static G4UnitsTable& GetUnitsTable();

  static G4bool IsUnitDefined(const G4String&);
  static G4double GetValueOf(const G4String&);
  static G4String GetCategory(const G4String&);

 private:
  G4UnitDefinition(const G4UnitDefinition&);
  G4UnitDefinition& operator=(const G4UnitDefinition&);

 private:
  G4String Name;         // SI name
  G4String SymbolName;   // SI symbol
  G4double Value = 0.0;  // value in the internal system of units

  static G4ThreadLocal G4UnitsTable* pUnitsTable;  // table of Units
  static G4ThreadLocal G4bool unitsTableDestroyed;

  std::size_t CategoryIndex = 0;  // category index of this unit

#ifdef G4MULTITHREADED
  static G4UnitsTable* pUnitsTableShadow;  // shadow of table of Units

 public:
  inline static G4UnitsTable& GetUnitsTableShadow()
  {
    return *pUnitsTableShadow;
  }
#endif
};

// --------------------------------------------------------------------

using G4UnitsContainer = std::vector<G4UnitDefinition*>;

class G4UnitsCategory
{
 public:
  explicit G4UnitsCategory(const G4String& name);
  ~G4UnitsCategory();
  G4bool operator==(const G4UnitsCategory&) const;
  G4bool operator!=(const G4UnitsCategory&) const;

  inline const G4String& GetName() const;
  inline G4UnitsContainer& GetUnitsList();
  inline G4int GetNameMxLen() const;
  inline G4int GetSymbMxLen() const;
  inline void UpdateNameMxLen(G4int len);
  inline void UpdateSymbMxLen(G4int len);
  void PrintCategory();

 private:
  G4UnitsCategory(const G4UnitsCategory&);
  G4UnitsCategory& operator=(const G4UnitsCategory&);

 private:
  G4String Name;               // dimensional family: Length,Volume,Energy
  G4UnitsContainer UnitsList;  // List of units in this family
  G4int NameMxLen = 0;         // max length of the units name
  G4int SymbMxLen = 0;         // max length of the units symbol
};

// --------------------------------------------------------------------

class G4BestUnit
{
 public:
  G4BestUnit(G4double internalValue, const G4String& category);
  G4BestUnit(const G4ThreeVector& internalValue, const G4String& category);
  // These constructors convert a physical quantity from its internalValue
  // into the most appropriate unit of the same category.
  // In practice it builds an object VU = (newValue, newUnit)

  ~G4BestUnit() = default;

  inline G4double* GetValue();
  inline const G4String& GetCategory() const;
  inline std::size_t GetIndexOfCategory() const;
  operator G4String() const;  // Conversion to best string.

  friend std::ostream& operator<<(std::ostream&, const G4BestUnit& VU);
  // Default format to print the objet VU above.

 private:
  G4double Value[3];   // value in the internal system of units
  G4int nbOfVals = 0;  // G4double=1; G4ThreeVector=3
  G4String Category;   // dimensional family: Length,Volume,Energy ...
  std::size_t IndexOfCategory = 0;  // position of Category in UnitsTable
};

#include "G4UnitsTable.icc"

#endif
