// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsTable.hh,v 1.5 1999-11-17 18:43:28 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -----------------------------------------------------------------
//	GEANT 4 class header file 
//
//      ------------------- class G4UnitsTable -----------------
//
// 17-05-98: first version, M.Maire
// 13-10-98: Units and symbols printed in fxed length

// class description
//
// This class maintains a table of Units.
// A Unit has a name, a symbol, a value and belong to a category (i.e. its
// dimensional definition): Length, Time, Energy, etc...
// The Units are grouped by category. The TableOfUnits is a list of categories.
// The class BestUnit allow to convert automaticaly a physical quantity
// from its internal value into the most appropriate Unit of the same category.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4UnitsTable_HH
#define G4UnitsTable_HH

#include "globals.hh"
#include "g4rw/tpordvec.h"

class G4UnitsCategory;
typedef G4RWTPtrOrderedVector<G4UnitsCategory> G4UnitsTable;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4UnitDefinition
{
public:  // with description

    G4UnitDefinition(G4String name, G4String symbol,G4String category,
            G4double value);
	    
public:  // without description
	    
   ~G4UnitDefinition();
    G4int operator==(const G4UnitDefinition&) const;
    G4int operator!=(const G4UnitDefinition&) const;
    
private:

    G4UnitDefinition(G4UnitDefinition&);
    G4UnitDefinition& operator=(const G4UnitDefinition&);
   
public:  // with description

    G4String      GetName()       {return Name;};
    G4String      GetSymbol()     {return SymbolName;};
    G4double      GetValue()      {return Value;};
    
    void          PrintDefinition();
    
    static void BuildUnitsTable();    
    static void PrintUnitsTable();
    
    static  
    G4UnitsTable&   GetUnitsTable() {return theUnitsTable;};
        
    static G4double GetValueOf (G4String);
    static G4String GetCategory(G4String);

private:

    G4String Name;            // SI name
    G4String SymbolName;      // SI symbol
    G4double Value;           // value in the internal system of units
    
    static 
    G4UnitsTable theUnitsTable;   // table of Units
    
    
    size_t CategoryIndex;         // category index of this unit
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4RWTPtrOrderedVector<G4UnitDefinition> G4UnitsContainer;

class G4UnitsCategory
{
public:  //without description

    G4UnitsCategory(G4String name);
   ~G4UnitsCategory();
    G4int operator==(const G4UnitsCategory&) const;
    G4int operator!=(const G4UnitsCategory&) const;
    
private:

    G4UnitsCategory(G4UnitsCategory&);
    G4UnitsCategory& operator=(const G4UnitsCategory&);
   
public:  //without description

    G4String          GetName()      {return Name;};
    G4UnitsContainer& GetUnitsList() {return UnitsList;};
    G4int             GetNameMxLen() {return NameMxLen;};
    G4int             GetSymbMxLen() {return SymbMxLen;};
    void  UpdateNameMxLen(G4int len) {if (NameMxLen<len) NameMxLen=len;};
    void  UpdateSymbMxLen(G4int len) {if (SymbMxLen<len) SymbMxLen=len;};
    void PrintCategory();

private:

    G4String          Name;        // dimensional family: Length,Volume,Energy ...
    G4UnitsContainer  UnitsList;   // List of units in this family
    G4int             NameMxLen;   // max length of the units name
    G4int             SymbMxLen;   // max length of the units symbol
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4BestUnit
{
public:  //with description

    G4BestUnit(G4double internalValue, G4String category);
    //This constructor converts a physical quantity from its internalValue
    //into the most appropriate unit of the same category.
    //In practice it builds an object VU = (newValue, newUnit)
    
   ~G4BestUnit();
   
public:  //without description

    G4double  GetValue()           {return Value;};
    G4String  GetCategory()        {return Category;};
    size_t    GetIndexOfCategory() {return IndexOfCategory;};
    
public:  //with description 
   
    friend
    ostream&  operator<<(ostream&,G4BestUnit VU);
    //default format to print the objet VU above

private:

    G4double   Value;        // value in the internal system of units
    G4String   Category;     // dimensional family: Length,Volume,Energy ...
    size_t IndexOfCategory;  // position of Category in UnitsTable
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
