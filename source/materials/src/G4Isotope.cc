//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Isotope.cc,v 1.11 2001-11-29 15:19:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 14.09.01: fCountUse: nb of elements which use this isotope 
// 13.09.01: suppression of the data member fIndexInTable
// 17.07.01: migration to STL. M. Verderi.
// 03.05.01: flux.precision(prec) at begin/end of operator<<
// 29.01.97: Forbidden to create Isotope with Z<1 or N<Z, M.Maire
// 26.06.96: Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Isotope.hh"
#include "g4std/iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IsotopeTable G4Isotope::theIsotopeTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Create an isotope
//
G4Isotope::G4Isotope(const G4String& Name, G4int Z, G4int N, G4double A)
: fName(Name), fZ(Z), fN(N), fA(A), fCountUse(0)
{
  if (Z<1) G4Exception
    (" ERROR! It is not allowed to create an Isotope with Z < 1" );

  if (N<Z) G4Exception
    (" ERROR! Attempt to create an Isotope with N < Z !!!" );

  theIsotopeTable.push_back(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope::~G4Isotope()
{
  if (fCountUse != 0)
    G4cout << "--> warning from ~G4Isotope(): the isotope " << fName
           << " is still referenced by " << fCountUse << " G4Elements \n" 
	   << G4endl;
	     
  //remove this isotope from theIsotopeTable
  G4IsotopeTable::iterator iter  = theIsotopeTable.begin();
  while ((iter != theIsotopeTable.end())&&(*iter != this)) iter++;
  if (iter != theIsotopeTable.end()) theIsotopeTable.erase(iter);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope::G4Isotope(G4Isotope& right)
{
  *this = right;
  
  //insert this new isotope in table
  theIsotopeTable.push_back(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope & G4Isotope::operator=(const G4Isotope& right)
{
  if (this != &right)
  {
    fName = right.fName;
    fZ = right.fZ;
    fN = right.fZ;
    fA = right.fA;
  }
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Isotope::operator==(const G4Isotope &right) const
{
  return (this == (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Isotope::operator!=(const G4Isotope &right) const
{
  return (this != (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4std::ostream& operator<<(G4std::ostream& flux, G4Isotope* isotope)
{
#ifdef G4USE_STD_NAMESPACE
  G4std::ios::fmtflags mode = flux.flags();
  flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
#else
  long mode = flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
#endif
  long prec = flux.precision(3);
    
  flux
    << " Isotope: " << G4std::setw(5) << isotope->fName 
    << "   Z = " << G4std::setw(2)    << isotope->fZ 
    << "   N = " << G4std::setw(3)    << isotope->fN
    << "   A = " << G4std::setw(6) << G4std::setprecision(2) 
    << (isotope->fA)/(g/mole) << " g/mole";

  flux.precision(prec);       
  flux.setf(mode,G4std::ios::floatfield);       
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4std::ostream& operator<<(G4std::ostream& flux, G4Isotope& isotope)
{
  flux << &isotope;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     
G4std::ostream& operator<<(G4std::ostream& flux, G4IsotopeTable IsotopeTable)
{
 //Dump info for all known isotopes
   flux 
     << "\n***** Table : Nb of isotopes = " << IsotopeTable.size() 
     << " *****\n" << G4endl;
        
   for (size_t i=0; i<IsotopeTable.size(); i++)
     flux << IsotopeTable[i] << G4endl;

   return flux;
}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4IsotopeTable* G4Isotope::GetIsotopeTable()
{
  return &theIsotopeTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t G4Isotope::GetNumberOfIsotopes()
{
  return theIsotopeTable.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope* G4Isotope::GetIsotope(G4String isotopeName)
{  
  // search the isotope by its name 
  for (size_t J=0 ; J<theIsotopeTable.size() ; J++)
   {
    if (theIsotopeTable[J]->GetName() == isotopeName)
      return theIsotopeTable[J];
   }
   
  // the isotope does not exist in the table
  return 0;          
}
