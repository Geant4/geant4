// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Isotope.cc,v 1.1 1999-01-07 16:09:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//
//      ---------- class G4Isotope ---------
//
//           Torre Wenaus, November 1995
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 26-06-96, Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban
// 29-01-97, Forbidden to create Isotope with Z<1 or N<Z, M.Maire

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4Isotope.hh"
#include <iomanip.h>

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

G4Isotope::G4Isotope(G4Isotope &right)
{
    *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Isotope & G4Isotope::operator=(const G4Isotope &right)
{
  if (this != &right)
  {
    fName = right.fName;
    fZ = right.fZ;
    fN = right.fZ;
    fA = right.fA;
    theIsotopeTable = right.theIsotopeTable;
    fIndexInTable = right.fIndexInTable;
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

ostream& operator<<(ostream& flux, G4Isotope* isotope)
{
  long mode = flux.setf(ios::fixed,ios::floatfield);
  
  flux
    << " Isotope: " << setw(5) << isotope->fName 
    << "   Z = " << setw(2) <<  isotope->fZ 
    << "   N = " << setw(3) <<  isotope->fN
    << "   A = " << setw(6) << setprecision(2) 
    << (isotope->fA)/(g/mole) << " g/mole";
    
  flux.setf(mode,ios::floatfield);       
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

 ostream& operator<<(ostream& flux, G4Isotope& isotope)
{
  flux << &isotope;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
     
ostream& operator<<(ostream& flux, G4IsotopeTable IsotopeTable)
{
 //Dump info for all known isotopes
   flux 
     << "\n***** Table : Nb of isotopes = " << IsotopeTable.length() 
     << " *****\n" << endl;
        
   for (G4int i=0; i<IsotopeTable.length(); i++) flux << IsotopeTable[i] << endl;

   return flux;
}      
