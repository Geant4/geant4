// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Isotope.cc,v 1.3 2001-03-12 17:48:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 29-01-97: Forbidden to create Isotope with Z<1 or N<Z, M.Maire
// 26-06-96: Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4Isotope.hh"
#include "g4std/iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4IsotopeTable G4Isotope::theIsotopeTable;

// Create an isotope
G4Isotope::G4Isotope(const G4String& Name, G4int Z, G4int N, G4double A)
: fName(Name), fZ(Z), fN(N), fA(A)
{
    if (Z<1) G4Exception
      (" ERROR! It is not allowed to create an Isotope with Z < 1" );

    if (N<Z) G4Exception
      (" ERROR! Attempt to create an Isotope with N < Z !!!" );

    theIsotopeTable.insert(this);
    fIndexInTable = theIsotopeTable.index(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Isotope::~G4Isotope() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Isotope::G4Isotope(G4Isotope& right)
{
  *this = right;
  
  //insert this new isotope in table
  theIsotopeTable.insert(this);
  fIndexInTable = theIsotopeTable.index(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Isotope & G4Isotope::operator=(const G4Isotope& right)
{
  if (this != &right)
  {
    fName = right.fName;
    fZ = right.fZ;
    fN = right.fZ;
    fA = right.fA;
    fIndexInTable = theIsotopeTable.index(this);     
  }
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4Isotope::operator==(const G4Isotope &right) const
{
  return (this == (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4Isotope::operator!=(const G4Isotope &right) const
{
  return (this != (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4std::ostream& operator<<(G4std::ostream& flux, G4Isotope* isotope)
{
  long mode = flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
  
  flux
    << " Isotope: " << G4std::setw(5) << isotope->fName 
    << "   Z = " << G4std::setw(2) <<  isotope->fZ 
    << "   N = " << G4std::setw(3) <<  isotope->fN
    << "   A = " << G4std::setw(6) << G4std::setprecision(2) 
    << (isotope->fA)/(g/mole) << " g/mole";
    
  flux.setf(mode,G4std::ios::floatfield);       
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

 G4std::ostream& operator<<(G4std::ostream& flux, G4Isotope& isotope)
{
  flux << &isotope;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
     
G4std::ostream& operator<<(G4std::ostream& flux, G4IsotopeTable IsotopeTable)
{
 //Dump info for all known isotopes
   flux 
     << "\n***** Table : Nb of isotopes = " << IsotopeTable.length() 
     << " *****\n" << G4endl;
        
   for (G4int i=0; i<IsotopeTable.length(); i++) flux << IsotopeTable[i] << G4endl;

   return flux;
}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
