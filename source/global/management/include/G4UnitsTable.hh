
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsTable.hh,v 1.3 1999-11-11 10:47:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -----------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//
//      ------------------- class G4UnitsTable -----------------
//
// 17-05-98: first version, M.Maire
// 13-10-98: Units and symbols printed in fxed length

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
public:
    G4UnitDefinition(G4String name, G4String symbol,G4String category,
            G4double value);
   ~G4UnitDefinition();
    G4int operator==(const G4UnitDefinition&) const;
    G4int operator!=(const G4UnitDefinition&) const;
    
private:
    G4UnitDefinition(G4UnitDefinition&);
    G4UnitDefinition& operator=(const G4UnitDefinition&);
   
public:
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
public:
    G4UnitsCategory(G4String name);
   ~G4UnitsCategory();
    G4int operator==(const G4UnitsCategory&) const;
    G4int operator!=(const G4UnitsCategory&) const;
    
private:
    G4UnitsCategory(G4UnitsCategory&);
    G4UnitsCategory& operator=(const G4UnitsCategory&);
   
public:
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
public:
    G4BestUnit(G4double,G4String);
   ~G4BestUnit();
   
public:
    G4double  GetValue()           {return Value;};
    G4String  GetCategory()        {return Category;};
    size_t    GetIndexOfCategory() {return IndexOfCategory;};
    
    friend
    ostream&  operator<<(ostream&,G4BestUnit);
    

private:
    G4double   Value;        // value in the internal system of units
    G4String   Category;     // dimensional family: Length,Volume,Energy ...
    size_t IndexOfCategory;  // position of Category in UnitsTable
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
