// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Element.hh,v 1.5 1999-12-16 18:11:08 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class description
//
// An element is a chemical element either directly defined in terms of
// its charactaristics: its name, symbol,
//                      Z (effective atomic number)
//                      N (effective number of nucleons)
//                      A (effective mass of a mole)
// or in terms of a collection of constituent isotopes with specified 
// relative abundance (i.e. fraction of nb of atomes per volume).
//
// Quantities, with physical meaning or not, which are constant in a given 
// element are computed and stored here as Derived data members.
//
// The class contains as a private static member the table of defined
// elements (an ordered vector of elements).
//
// Elements can be assembled singly or in mixtures into materials used
// in volume definitions via the G4Material class.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 09-07-96, new data members added by L.Urban
// 17-01-97, aesthetic rearrangement, M.Maire
// 20-01-97, Tsai formula for the rad length, M.Maire
// 21-01-97, remove mixture flag, M.Maire
// 24-01-97, new data member: fTaul
//           new method: ComputeIonisationPara, M.Maire
// 20-03-97, corrected initialization of pointers, M.Maire
// 27-06-97, new function GetIsotope(int), M.Maire
// 24-02-98, fWeightVector becomes fRelativeAbundanceVector
// 27-04-98, atomic shell stuff, V. Grichine
// 09-07-98, Ionisation parameters removed from the class, M.Maire
// 04-08-98, new method GetElement(elementName), M.Maire
// 16-11-98, Subshell -> Shell, mma

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4ELEMENT_HH
#define G4ELEMENT_HH

#include "G4ios.hh"
#include "g4rw/tpvector.h"
#include "g4rw/tpordvec.h"
#include "globals.hh"
#include "G4Isotope.hh"
#include "G4AtomicShells.hh"
#include "G4IonisParamElm.hh"

typedef G4RWTPtrVector<G4Isotope> G4IsotopeVector;

class G4Element;              //forward declaration
typedef G4RWTPtrOrderedVector<G4Element> G4ElementTable;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4Element
{
public:  // with description

    //
    // Constructor to Build an element directly; no reference to isotopes
    //
    G4Element(const G4String& name,		//its name
              const G4String& symbol,		//its symbol
                    G4double  Zeff,		//atomic number
                    G4double  Aeff);		//mass of mole
                    
    //
    // Constructor to Build an element from isotopes via AddIsotope
    //
    G4Element(const G4String& name,		//its name
              const G4String& symbol,		//its symbol
              G4int nbIsotopes);		//nb of isotopes

    //
    // Add an isotope to the element
    // 
    void AddIsotope(G4Isotope* isotope,			//isotope 
                    G4double   RelativeAbundance);	//fraction of nb of 
                    					//atomes per volume
                    					
   ~G4Element();
   
    //
    // retrieval methods
    //
    G4String GetName()   const {return fName;};
    G4String GetSymbol() const {return fSymbol;};
    G4double GetZ()      const {return fZeff;};     //atomic number
    G4double GetN()      const {return fNeff;};     //number of nucleons
    G4double GetA()      const {return fAeff;};     //mass of a mole
    
    //the number of atomic shells in this element:
    G4int    GetNbOfAtomicShells() const {return fNbOfAtomicShells;};
    //the binding energy of the shell :
    G4double GetAtomicShell(G4int) const;
    
    //number of isotopes constituing this element:
    size_t  GetNumberOfIsotopes() const {return fNumberOfIsotopes;};
    
    //vector of pointers to isotopes constituing this element:
    G4IsotopeVector* GetIsotopeVector()    const {return theIsotopeVector;};
    
    //vector of relative abundance of each isotope:
    G4double* GetRelativeAbundanceVector() const {return fRelativeAbundanceVector;};
    
    const G4Isotope* GetIsotope(G4int iso) const {return (*theIsotopeVector)[iso];};

    //the (static) Table of Elements:
    static const  G4ElementTable* GetElementTable() {return &theElementTable;};
    static size_t GetNumberOfElements() {return theElementTable.length();};
    //the index of this element in the Table:
    size_t GetIndex() const             {return fIndexInTable;};
    
    //return  pointer to an element, given its name:
    static  G4Element* GetElement(G4String name);
    
    //Coulomb correction factor:
    G4double GetfCoulomb() const {return fCoulomb;};
    
    //Tsai formula for the radiation length:
    G4double GetfRadTsai() const {return fRadTsai;};
    
    //pointer to ionisation parameters:
    G4IonisParamElm* GetIonisation() const {return fIonisation;};
    
    //
    // printing methods
    //    
    friend G4std::ostream& operator<<(G4std::ostream&, G4Element*);    
    friend G4std::ostream& operator<<(G4std::ostream&, G4Element&);    
    friend G4std::ostream& operator<<(G4std::ostream&, G4ElementTable);

public:  // without description
    
    G4int operator==(const G4Element&) const;
    G4int operator!=(const G4Element&) const;
     
private:

    G4Element(G4Element &right);
    const G4Element & operator=(const G4Element &right);

private:

    void InitializePointers();
    void ComputeDerivedQuantities();
    void ComputeCoulombFactor();
    void ComputeLradTsaiFactor();

private:

  //
  // Basic data members (which define an Element)
  //
    G4String fName;              // name
    G4String fSymbol;            // symbol
    G4double fZeff;              // Effective atomic number
    G4double fNeff;              // Effective number of nucleons
    G4double fAeff;              // Effective mass of a mole
    
    G4int fNbOfAtomicShells;     // number  of atomic shells
    G4double* fAtomicShells ;    // Pointer to atomic shell binding energies
    
    // Isotope vector contains constituent isotopes of the element   
    size_t fNumberOfIsotopes;    // Number of isotopes added to the element
    G4IsotopeVector* theIsotopeVector;
    G4double* fRelativeAbundanceVector;     // Fraction of nb of atomes per volume
                                            // for each constituent

    // Set up the static Table of Elements
    static G4ElementTable theElementTable;
    size_t fIndexInTable;          // Position of the element in the table

  //
  // Derived data members (computed from the basic data members)
  //
    G4double fCoulomb;             // Coulomb correction factor
    G4double fRadTsai;             // Tsai formula for the radiation length
    G4IonisParamElm* fIonisation;  // Pointer to ionisation parameters
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4Element* G4Element::GetElement(G4String elementName)
{  
  // search the element by its name 
  for (G4int J=0 ; J<theElementTable.length() ; J++)
   {
    if(theElementTable[J]->GetName() == elementName)
      return theElementTable[J];
   }
   
  G4cerr << "  Warning from GetElement(name). The element: " << elementName
         << "  does not exist in the ElementTable.  Return NULL pointer \n";
  return NULL;   
}

#endif
