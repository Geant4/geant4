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

//---------------------------------------------------------------------------
//
// ClassName:   G4Element
//
// Description: Contains element properties
//
// Class description:
//
// An element is a chemical element either directly defined in terms of
// its characteristics: its name, symbol,
//                      Z (effective atomic number)
//                      N (effective number of nucleons)
//                      A (effective mass of a mole)
// or in terms of a collection of constituent isotopes with specified 
// relative abundance (i.e. fraction of nb of atoms per volume).
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
// It is strongly recommended do not delete G4Element instance in the
// user code. All G4Elements will be automatically deleted at the end 
// of Geant4 session

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// 30-03-01, suppression of the warning message in GetElement
// 17-07-01, migration to STL, M. Verderi
// 13-09-01, stl migration. Suppression of the data member fIndexInTable
// 14-09-01, fCountUse: nb of materials which use this element
// 26-02-02, fIndexInTable renewed 
// 01-04-05, new data member fIndexZ to count the number of elements with same Z
// 17-10-06: Add Get/Set fNaturalAbundance (V.Ivanchenko)
// 17.09.09, add fNbOfShellElectrons and methods (V. Grichine)
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4ELEMENT_HH
#define G4ELEMENT_HH 1

#include "globals.hh"
#include <vector>
#include "G4ios.hh"
#include "G4Isotope.hh"
#include "G4IonisParamElm.hh"
#include "G4IsotopeVector.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
            G4int nbIsotopes);			//nb of isotopes

  //
  // Add an isotope to the element
  // 
  void AddIsotope(G4Isotope* isotope,			//isotope 
                  G4double   RelativeAbundance);	//fraction of nb of 
                  					//atomes per volume
  virtual ~G4Element();
  
  //
  // Retrieval methods
  //
  inline const G4String& GetName()   const {return fName;}
  inline const G4String& GetSymbol() const {return fSymbol;}

  // Atomic number
  inline G4double GetZ()             const {return fZeff;}    
  inline G4int GetZasInt()           const {return fZ;}    

  // Atomic weight in atomic units
  inline G4double GetN()             const {return fNeff;}     
  inline G4double GetAtomicMassAmu() const {return fNeff;}      

  // Mass of a mole in Geant4 units for atoms with atomic shell
  inline G4double GetA()             const {return fAeff;}    

  inline G4bool   GetNaturalAbundanceFlag() const;

  inline void     SetNaturalAbundanceFlag(G4bool);
  
  //the number of atomic shells in this element:
  //
  inline G4int GetNbOfAtomicShells() const {return fNbOfAtomicShells;}
  
  //the binding energy of the shell, ground shell index=0
  //
  G4double GetAtomicShell(G4int index) const;

  //the number of electrons at the shell, ground shell index=0
  //
  G4int GetNbOfShellElectrons(G4int index) const;
    
  //number of isotopes constituing this element:
  //
  inline size_t GetNumberOfIsotopes() const {return fNumberOfIsotopes;}
   
  //vector of pointers to isotopes constituing this element:
  //
  inline G4IsotopeVector* GetIsotopeVector() const {return theIsotopeVector;}
    
  //vector of relative abundance of each isotope:
  //
  inline G4double* GetRelativeAbundanceVector() const 
                   {return fRelativeAbundanceVector;}
    
  inline const G4Isotope* GetIsotope(G4int iso) const 
                   {return (*theIsotopeVector)[iso];}

  //the (static) Table of Elements:
  //
  static G4ElementTable* GetElementTable();
  
  static 
  size_t GetNumberOfElements();
  
  //the index of this element in the Table:
  //
  inline size_t GetIndex() const {return fIndexInTable;}
    
  //return pointer to an element, given its name:
  //
  static G4Element* GetElement(const G4String& name, G4bool warning = true);

  //Coulomb correction factor:
  //
  inline G4double GetfCoulomb() const {return fCoulomb;}
   
  //Tsai formula for the radiation length:
  //
  inline G4double GetfRadTsai() const {return fRadTsai;}
    
  //pointer to ionisation parameters:
  //
  inline G4IonisParamElm* GetIonisation() const {return fIonisation;}
    
  // printing methods
  //    
  friend std::ostream& operator<<(std::ostream&, const G4Element*);
  friend std::ostream& operator<<(std::ostream&, const G4Element&);
  friend std::ostream& operator<<(std::ostream&, const G4ElementTable&);
  friend std::ostream& operator<<(std::ostream&, const G4ElementVector&);

public:  // without description

  G4Element(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  inline void SetName(const G4String& name)  {fName=name;}

  G4Element(G4Element&) = delete;
  const G4Element & operator=(const G4Element&) = delete;
  G4bool operator==(const G4Element&) const = delete;
  G4bool operator!=(const G4Element&) const = delete;

private:

  void InitializePointers();
  void ComputeDerivedQuantities();
  void ComputeCoulombFactor();
  void ComputeLradTsaiFactor();
  void AddNaturalIsotopes();

  //
  // Basic data members (which define an Element)
  //
  G4String fName;              // name
  G4String fSymbol;            // symbol
  G4double fZeff;              // Effective atomic number
  G4double fNeff;              // Effective number of nucleons
  G4double fAeff;              // Effective mass of a mole
  G4int    fZ;
    
  G4int fNbOfAtomicShells;     // number  of atomic shells
  G4double* fAtomicShells ;    // Pointer to atomic shell binding energies
  G4int* fNbOfShellElectrons;  // Pointer to the number of subshell electrons
    
  // Isotope vector contains constituent isotopes of the element   
  G4int fNumberOfIsotopes;     // Number of isotopes added to the element
  G4IsotopeVector* theIsotopeVector;
  G4double* fRelativeAbundanceVector;     // Fraction nb of atomes per volume
                                          // for each constituent

  // Set up the static Table of Elements
  static G4ElementTable theElementTable;
  size_t fIndexInTable;
  G4bool fNaturalAbundance;

  //
  // Derived data members (computed from the basic data members)
  //
  G4double fCoulomb;             // Coulomb correction factor
  G4double fRadTsai;             // Tsai formula for the radiation length
  G4IonisParamElm* fIonisation;  // Pointer to ionisation parameters
};

inline G4bool G4Element::GetNaturalAbundanceFlag() const
{
  return fNaturalAbundance;
}

inline void G4Element::SetNaturalAbundanceFlag(G4bool val) 
{
  fNaturalAbundance = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
