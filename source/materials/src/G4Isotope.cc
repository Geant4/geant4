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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 26.06.96: Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban
// 29.01.97: Forbidden to create Isotope with Z<1 or N<Z, M.Maire
// 03.05.01: flux.precision(prec) at begin/end of operator<<
// 17.07.01: migration to STL. M. Verderi.
// 13.09.01: suppression of the data member fIndexInTable
// 14.09.01: fCountUse: nb of elements which use this isotope 
// 26.02.02: fIndexInTable renewed
// 17.10.06: if fA is not defined in the constructor, it is computed from
//           NistManager v.Ivanchenko
// 25.10.11: new scheme for G4Exception  (mma)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iomanip>

#include "G4Isotope.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4IsotopeTable G4Isotope::theIsotopeTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Create an isotope
//
G4Isotope::G4Isotope(const G4String& Name,G4int Z,G4int N,G4double A,G4int il)
  : fName(Name), fZ(Z), fN(N), fA(A), fm(il)
{
  if (Z<1) { 
    G4ExceptionDescription ed;
    ed << "Wrong Isotope " << Name << " Z= " << Z << G4endl;
    G4Exception ("G4Isotope::G4Isotope()", "mat001", FatalException, ed);
  }
  if (N<Z) {
    G4ExceptionDescription ed;
    ed << "Wrong Isotope " << Name << " Z= " << Z << " > N= " << N << G4endl;
    G4Exception ("G4Isotope::G4Isotope()", "mat002", FatalException, ed);
  }
  if (A<=0.0) {
    fA = (G4NistManager::Instance()->GetAtomicMass(Z,N))
      *CLHEP::g/(CLHEP::mole*CLHEP::amu_c2);  
  }
  theIsotopeTable.push_back(this);
  fIndexInTable = theIsotopeTable.size() - 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Isotope::G4Isotope(__void__&)
  : fZ(0), fN(0), fA(0), fm(0), fIndexInTable(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope::~G4Isotope()
{
  theIsotopeTable[fIndexInTable] = nullptr;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope::G4Isotope(G4Isotope& right)
{
  *this = right;
  
  //insert this new isotope in table
  theIsotopeTable.push_back(this);
  fIndexInTable = theIsotopeTable.size() - 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Isotope & G4Isotope::operator=(const G4Isotope& right)
{
  if (this != &right)
  {
    fName = right.fName;
    fZ = right.fZ;
    fN = right.fN;
    fA = right.fA;
    fm = right.fm;    
  }
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4Isotope::operator==(const G4Isotope &right) const
{
  return (this == (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4Isotope::operator!=(const G4Isotope &right) const
{
  return (this != (G4Isotope *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream& operator<<(std::ostream& flux, const G4Isotope* isotope)
{
  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  G4long prec = flux.precision(3);
    
  flux
    << " Isotope: " << std::setw(5) << isotope->fName 
    << "   Z = " << std::setw(2)    << isotope->fZ 
    << "   N = " << std::setw(3)    << isotope->fN
    << "   A = " << std::setw(6) << std::setprecision(2) 
    << (isotope->fA)/(g/mole) << " g/mole";

  flux.precision(prec);       
  flux.setf(mode,std::ios::floatfield);       
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream& operator<<(std::ostream& flux, const G4Isotope& isotope)
{
  flux << &isotope;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     
std::ostream& operator<<(std::ostream& flux, G4IsotopeTable IsotopeTable)
{
 //Dump info for all known isotopes
   flux 
     << "\n***** Table : Nb of isotopes = " << IsotopeTable.size() 
     << " *****\n" << G4endl;

   for(auto& i : IsotopeTable)
   {
     flux << i << G4endl;
   }

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

G4Isotope* G4Isotope::GetIsotope(const G4String& isotopeName, G4bool warning)
{  
  // search the isotope by its name
  for(auto& J : theIsotopeTable)
  {
    if(J->GetName() == isotopeName)
    {
      return J;
    }
   }
   
  // the isotope does not exist in the table
  if (warning) {
    G4cout << "\n---> warning from G4Isotope::GetIsotope(). The isotope: "
	   << isotopeName << " does not exist in the table. Return NULL pointer."
	   << G4endl;
  }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
