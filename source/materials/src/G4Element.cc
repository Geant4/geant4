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
// $Id: G4Element.cc,v 1.14 2001-11-29 15:19:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 26-06-96: Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban
// 09-07-96: new data members added by L.Urban
// 17-01-97: aesthetic rearrangement, M.Maire
// 20-01-97: Compute Tsai's formula for the rad length, M.Maire
// 21-01-97: remove mixture flag, M.Maire
// 24-01-97: ComputeIonisationParameters(). 
//           new data member: fTaul, M.Maire
// 29-01-97: Forbidden to create Element with Z<1 or N<Z, M.Maire
// 20-03-97: corrected initialization of pointers, M.Maire
// 28-04-98: atomic subshell binding energies stuff, V. Grichine  
// 09-07-98: Ionisation parameters removed from the class, M.Maire
// 16-11-98: name Subshell -> Shell; GetBindingEnergy() (mma)
// 09-03-01: assignement operator revised (mma)
// 02-05-01: check identical Z in AddIsotope (marc)
// 03-05-01; flux.precision(prec) at begin/end of operator<<
// 13-09-01: suppression of the data member fIndexInTable
// 14-09-01: fCountUse: nb of materials which use this element
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Element.hh"
#include "g4std/iomanip"

G4ElementTable G4Element::theElementTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to Generate an element from scratch

G4Element::G4Element(const G4String& name, const G4String& symbol,
                     G4double zeff, G4double aeff)
:fName(name),fSymbol(symbol)		     
{
    if (zeff<1.) G4Exception (" ERROR from G4Element::G4Element !"
       " It is not allowed to create an Element with Z < 1" );

    if (aeff/(g/mole)<zeff) G4Exception (" ERROR from G4Element::G4Element !"
       " Attempt to create an Element with N < Z !!!" );

    if ((zeff-G4int(zeff)) > perMillion)
      G4cerr << name << " : WARNING from G4Element::G4Element !"  
         " Trying to define an element as a mixture directly via effective Z."
         << G4endl;

    InitializePointers();

    fZeff   = zeff;
    fNeff   = aeff/(g/mole);
    fAeff   = aeff;

    fNumberOfIsotopes = 0;
   
    fNbOfAtomicShells = G4AtomicShells::GetNumberOfShells((G4int)fZeff);
    fAtomicShells     = new G4double[fNbOfAtomicShells];
    for (G4int i=0;i<fNbOfAtomicShells;i++)
       fAtomicShells[i] = G4AtomicShells::GetBindingEnergy((G4int)fZeff,i);
           
    ComputeDerivedQuantities();

    // Store in ElementTable
    theElementTable.push_back(this);
    
    fCountUse = 0;           //nb of materials which use this element
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to Generate element from a List of 'nIsotopes' isotopes, added
// via AddIsotope

G4Element::G4Element(const G4String& name, const G4String& symbol, G4int nIsotopes)
:fName(name),fSymbol(symbol)
{
    InitializePointers();

    size_t n = size_t(nIsotopes);

    fNumberOfIsotopes = 0;

    theIsotopeVector         = new G4IsotopeVector(n,0);
    fRelativeAbundanceVector = new G4double[nIsotopes];
    
    fCountUse = 0;           //nb of materials which use this element    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add an isotope to the element

void G4Element::AddIsotope(G4Isotope* isotope, G4double abundance)
{
    if (theIsotopeVector == 0)
       G4Exception ("ERROR from G4Element::AddIsotope!"
       " Trying to add an Isotope before contructing the element.");

    // filling ...
    if ( fNumberOfIsotopes < theIsotopeVector->size() ) {
       // check same Z
       if (fNumberOfIsotopes==0) fZeff = G4double(isotope->GetZ());
       else if (G4double(isotope->GetZ()) != fZeff) 
          G4Exception ("ERROR from G4Element::AddIsotope!"
	   " Try to add isotopes with different Z");
       //Z ok   
       fRelativeAbundanceVector[fNumberOfIsotopes] = abundance;
       (*theIsotopeVector)[fNumberOfIsotopes] = isotope;
       ++fNumberOfIsotopes;
       isotope->increaseCountUse();
      } 
    else G4Exception ("ERROR from G4Element::AddIsotope!"  
       " Attempt to add more than the declared number of isotopes.");

    // filled.
    if ( fNumberOfIsotopes == theIsotopeVector->size() ) {
      // Compute Neff, Aeff
      G4double wtSum=0.0;

      fNeff = fAeff = 0.0;
      for (size_t i=0;i<fNumberOfIsotopes;i++) {
	fNeff +=  fRelativeAbundanceVector[i]*(*theIsotopeVector)[i]->GetN();
	fAeff +=  fRelativeAbundanceVector[i]*(*theIsotopeVector)[i]->GetA();
        wtSum +=  fRelativeAbundanceVector[i];
      }
      fNeff /=  wtSum;
      fAeff /=  wtSum;
      
      fNbOfAtomicShells = G4AtomicShells::GetNumberOfShells((G4int)fZeff);
      fAtomicShells     = new G4double[fNbOfAtomicShells];
      for (G4int j=0;j<fNbOfAtomicShells;j++)
         fAtomicShells[j] = G4AtomicShells::GetBindingEnergy((G4int)fZeff,j);
         
      ComputeDerivedQuantities();

      // Store in table
      theElementTable.push_back(this);
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Element::InitializePointers()
{
    theIsotopeVector = 0;
    fRelativeAbundanceVector = 0;
    fAtomicShells = 0;
    fIonisation = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Element::~G4Element()
{
  if (fCountUse != 0)
    G4cout << "--> warning from ~G4Element(): the element " << fName
           << " is still referenced by " << fCountUse << " G4Materials \n" 
	   << G4endl;
	   
  if (theIsotopeVector)
    { for (size_t i=0; i<fNumberOfIsotopes; i++)
                          (*theIsotopeVector)[i]->decreaseCountUse();         
      delete theIsotopeVector;
    } 
  if (fRelativeAbundanceVector) delete [] fRelativeAbundanceVector;
  if (fAtomicShells)            delete [] fAtomicShells;
  if (fIonisation)              delete    fIonisation;
  
  //remove this element from theElementTable
  G4ElementTable::iterator iter  = theElementTable.begin();
  while ((iter != theElementTable.end())&&(*iter != this)) iter++;
  if (iter != theElementTable.end()) theElementTable.erase(iter);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Element::ComputeDerivedQuantities()
{
  // some basic functions of the atomic number

  // Radiation Length
     ComputeCoulombFactor();
     ComputeLradTsaiFactor(); 

  // parameters for energy loss by ionisation 
     if (fIonisation) delete fIonisation;  
     fIonisation = new G4IonisParamElm(fZeff);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Element::ComputeCoulombFactor()
{
  //
  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const G4double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;

  G4double az2 = (fine_structure_const*fZeff)*(fine_structure_const*fZeff);
  G4double az4 = az2 * az2;

  fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Element::ComputeLradTsaiFactor()
{
  //
  //  Compute Tsai's Expression for the Radiation Length
  //  (Phys Rev. D50 3-1 (1994) page 1254)

  const G4double Lrad_light[]  = {5.31  , 4.79  , 4.74 ,  4.71} ;
  const G4double Lprad_light[] = {6.144 , 5.621 , 5.805 , 5.924} ;
  
  const G4double logZ3 = log(fZeff)/3.;

  G4double Lrad, Lprad;
  G4int iz = (int)(fZeff+0.5) - 1 ;
  if (iz <= 3) { Lrad = Lrad_light[iz] ;  Lprad = Lprad_light[iz] ; }
     else { Lrad = log(184.15) - logZ3 ; Lprad = log(1194.) - 2*logZ3 ; }

  fRadTsai = 4*alpha_rcl2*fZeff*(fZeff*(Lrad-fCoulomb) + Lprad); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4Element::GetAtomicShell(G4int i) const
{
  if (i<0 || i>=fNbOfAtomicShells)
      G4Exception("Invalid argument in G4Element::GetAtomicShell");
  return fAtomicShells[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ElementTable* G4Element::GetElementTable()
{
  return &theElementTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t G4Element::GetNumberOfElements()
{
  return theElementTable.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Element* G4Element::GetElement(G4String elementName)
{  
  // search the element by its name 
  for (size_t J=0 ; J<theElementTable.size() ; J++)
   {
    if (theElementTable[J]->GetName() == elementName)
      return theElementTable[J];
   }
   
  // the element does not exist in the table 
  return 0;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Element::G4Element(G4Element& right)
{
      InitializePointers();
      *this = right;

      // Store this new element in table and set the index
      theElementTable.push_back(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4Element& G4Element::operator=(const G4Element& right)
{
  if (this != &right)
    {
      fName                    = right.fName;
      fSymbol                  = right.fSymbol;
      fZeff                    = right.fZeff;
      fNeff                    = right.fNeff;
      fAeff                    = right.fAeff;
      
      if (fAtomicShells) delete [] fAtomicShells;      
      fNbOfAtomicShells        = right.fNbOfAtomicShells;
      fAtomicShells     = new G4double[fNbOfAtomicShells];
      for (G4int i=0;i<fNbOfAtomicShells;i++)      
         fAtomicShells[i]      = right.fAtomicShells[i];
	 
      if (theIsotopeVector) delete theIsotopeVector;
      if (fRelativeAbundanceVector) delete [] fRelativeAbundanceVector;
	      	 
      fNumberOfIsotopes        = right.fNumberOfIsotopes;
      if (fNumberOfIsotopes > 0)
        {
	 theIsotopeVector         = new G4IsotopeVector(fNumberOfIsotopes,0);
	 fRelativeAbundanceVector = new G4double[fNumberOfIsotopes];
	 for (size_t i=0;i<fNumberOfIsotopes;i++)
	    {
             (*theIsotopeVector)[i]      = (*right.theIsotopeVector)[i];
             fRelativeAbundanceVector[i] = right.fRelativeAbundanceVector[i];
	    }
	}   
      ComputeDerivedQuantities();
     } 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Element::operator==(const G4Element& right) const
{
  return (this == (G4Element*) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4Element::operator!=(const G4Element& right) const
{
  return (this != (G4Element*) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4std::ostream& operator<<(G4std::ostream& flux, G4Element* element)
{
#ifdef G4USE_STD_NAMESPACE
  G4std::ios::fmtflags mode = flux.flags();
  flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
#else 
  long mode = flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
#endif
  long prec = flux.precision(3);
  
  flux
    << " Element: " << G4std::setw(8) << element->fName << G4std::setw(3)
                    << element->fSymbol
    << "   Z = " << G4std::setw(4) << G4std::setprecision(1) <<  element->fZeff 
    << "   N = " << G4std::setw(5) << G4std::setprecision(1) <<  element->fNeff
    << "   A = " << G4std::setw(6) << G4std::setprecision(2)
                 << (element->fAeff)/(g/mole) << " g/mole";
   
  for (size_t i=0; i<element->fNumberOfIsotopes; i++)
  flux 
    << "\n   ---> " << (*(element->theIsotopeVector))[i] 
    << "   abundance: " << G4std::setw(6) << G4std::setprecision(2) 
    << (element->fRelativeAbundanceVector[i])/perCent << " %";
    
  flux.precision(prec);        
  flux.setf(mode,G4std::ios::floatfield);         
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 G4std::ostream& operator<<(G4std::ostream& flux, G4Element& element)
{
  flux << &element;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     
G4std::ostream& operator<<(G4std::ostream& flux, G4ElementTable ElementTable)
{
 //Dump info for all known elements
   flux << "\n***** Table : Nb of elements = " << ElementTable.size() 
        << " *****\n" << G4endl;
        
   for (size_t i=0; i<ElementTable.size(); i++) flux << ElementTable[i] 
						     << G4endl << G4endl;

   return flux;
}
      
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
