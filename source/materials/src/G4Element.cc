// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Element.cc,v 1.3 1999-12-15 14:50:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//      ------------ class G4Element ------------
//
//             Torre Wenaus, November 1995
//
//
// 26-06-96, Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban
// 09-07-96, new data members added by L.Urban
// 17-01-97, aesthetic rearrangement, M.Maire
// 20-01-97, Compute Tsai's formula for the rad length, M.Maire
// 21-01-97, remove mixture flag, M.Maire
// 24-01-97, ComputeIonisationParameters(). 
//           new data member: fTaul, M.Maire
// 29-01-97, Forbidden to create Element with Z<1 or N<Z, M.Maire
// 20-03-97, corrected initialization of pointers, M.Maire
// 28-04-98, atomic subshell binding energies stuff, V. Grichine  
// 09-07-98, Ionisation parameters removed from the class, M.Maire
// 16-11-98, name Subshell -> Shell; GetBindingEnergy(), mma

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#include "G4Element.hh"
#include "g4std/iomanip"

G4ElementTable G4Element::theElementTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Constructor to Generate an element from scratch

G4Element::G4Element(const G4String& name, const G4String& symbol,
                     G4double zeff, G4double aeff)
:fName(name),fSymbol(symbol)		     
{
    if (zeff<1.) G4Exception
      (" ERROR! It is not allowed to create an Element with Z < 1" );

    if (aeff/(g/mole)<zeff) G4Exception
      (" ERROR! Attempt to create an Element with N < Z !!!" );

    if ((zeff-G4int(zeff)) > perMillion)
      G4cerr << name <<
      " : WARNING ! Trying to define an element as a mixture directly via effective Z."
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

    // Store in table and set the index
    theElementTable.insert(this);
    fIndexInTable = theElementTable.index(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Constructor to Generate element from a List of 'nIsotopes' isotopes, added
// via AddIsotope

G4Element::G4Element(const G4String& name, const G4String& symbol, G4int nIsotopes)
:fName(name),fSymbol(symbol)
{
    InitializePointers();

    size_t n = size_t(nIsotopes);

    fNumberOfIsotopes = 0;

    theIsotopeVector         = new G4IsotopeVector(n);
    fRelativeAbundanceVector = new G4double[nIsotopes];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Add an isotope to the element

void G4Element::AddIsotope(G4Isotope* isotope, G4double abundance)
{
    if (theIsotopeVector == NULL)
       G4Exception("ERROR!!! - Trying to add an Isotope before contructing the element.");

    // filling ...
    if ( fNumberOfIsotopes < theIsotopeVector->length() ) {
       fRelativeAbundanceVector[fNumberOfIsotopes] = abundance;
       (*theIsotopeVector)(fNumberOfIsotopes) = isotope;
       fNumberOfIsotopes ++;
      } 
    else
       G4Exception
      ("ERROR!!! - Attempt to add more than the declared number of constituent isotopes.");

    // filled.
    if ( fNumberOfIsotopes == theIsotopeVector->length() ) {
      // Compute Zeff, Neff, Aeff
      G4int i;
      G4double wtSum=0.0;

      fZeff = G4double( (*theIsotopeVector)(0)->GetZ() );
      fNeff = fAeff = 0.0;
      for (i=0;i<fNumberOfIsotopes;i++) {
        fNeff +=  fRelativeAbundanceVector[i]*(*theIsotopeVector)(i)->GetN();
        fAeff +=  fRelativeAbundanceVector[i]*(*theIsotopeVector)(i)->GetA();
        wtSum +=  fRelativeAbundanceVector[i];
      }
      fNeff /=  wtSum;
      fAeff /=  wtSum;
      
      fNbOfAtomicShells = G4AtomicShells::GetNumberOfShells((G4int)fZeff);
      fAtomicShells     = new G4double[fNbOfAtomicShells];
      for (i=0;i<fNbOfAtomicShells;i++)
         fAtomicShells[i] = G4AtomicShells::GetBindingEnergy((G4int)fZeff,i);
         
      ComputeDerivedQuantities();

      // Store in table and set the index
      theElementTable.insert(this);
      fIndexInTable = theElementTable.index(this);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4Element::InitializePointers()
{
    theIsotopeVector = NULL;
    fRelativeAbundanceVector = NULL;
    fAtomicShells = NULL;
    fIonisation = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Element::~G4Element()
{
  if (theIsotopeVector)         delete    theIsotopeVector;
  if (fRelativeAbundanceVector) delete [] fRelativeAbundanceVector;
  if (fAtomicShells)            delete [] fAtomicShells;
  if (fIonisation)              delete    fIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4Element::G4Element(G4Element &right)
{
    *this = right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

const G4Element& G4Element::operator=(const G4Element& right)
{
  if (this != &right)
    {
      fName                    = right.fName;
      fSymbol                  = right.fSymbol;
      fZeff                    = right.fZeff;
      fNeff                    = right.fNeff;
      fAeff                    = right.fAeff;
      fNbOfAtomicShells        = right.fNbOfAtomicShells;
      fAtomicShells            = right.fAtomicShells;
      fNumberOfIsotopes        = right.fNumberOfIsotopes;
      theIsotopeVector         = right.theIsotopeVector;
      fRelativeAbundanceVector = right.fRelativeAbundanceVector;
      fIndexInTable            = right.fIndexInTable;
      fCoulomb                 = right.fCoulomb;
      fRadTsai                 = right.fRadTsai;
      fIonisation              = right.fIonisation;
     } 
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4Element::operator==(const G4Element &right) const
{
  return (this == (G4Element *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4int G4Element::operator!=(const G4Element &right) const
{
  return (this != (G4Element *) &right);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4Element::ComputeDerivedQuantities()
{
  // some basic functions of the atomic number

  // Radiation Length
     ComputeCoulombFactor();
     ComputeLradTsaiFactor(); 

  // parameters for energy loss by ionisation   
     fIonisation = new G4IonisParamElm(fZeff);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4Element::ComputeCoulombFactor()
{
  //
  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const G4double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;

  G4double az2 = (fine_structure_const*fZeff)*(fine_structure_const*fZeff);
  G4double az4 = az2 * az2;

  fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double G4Element::GetAtomicShell(G4int i) const
{
  if (i<0 || i>=fNbOfAtomicShells)
      G4Exception("Invalid argument in G4Element::GetAtomicShell");
  return fAtomicShells[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4std::ostream& operator<<(G4std::ostream& flux, G4Element* element)
{ 
  long mode = flux.setf(G4std::ios::fixed,G4std::ios::floatfield);
  
  flux
    << " Element: " << G4std::setw(8) << element->fName << G4std::setw(3) << element->fSymbol
    << "   Z = " << G4std::setw(4) << G4std::setprecision(1) <<  element->fZeff 
    << "   N = " << G4std::setw(5) << G4std::setprecision(1) <<  element->fNeff
    << "   A = " << G4std::setw(6) << G4std::setprecision(2) << (element->fAeff)/(g/mole) 
    << " g/mole";
   
  for (G4int i=0; i<element->fNumberOfIsotopes; i++)
  flux 
    << "\n   ---> " << (*(element->theIsotopeVector))[i] 
    << "   abundance: " << G4std::setw(6) << G4std::setprecision(2) 
    << (element->fRelativeAbundanceVector[i])/perCent << " %";
    
  flux.setf(mode,G4std::ios::floatfield);         
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

 G4std::ostream& operator<<(G4std::ostream& flux, G4Element& element)
{
  flux << &element;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
     
G4std::ostream& operator<<(G4std::ostream& flux, G4ElementTable ElementTable)
{
 //Dump info for all known elements
   flux << "\n***** Table : Nb of elements = " << ElementTable.length() 
        << " *****\n" << G4endl;
        
   for (G4int i=0; i<ElementTable.length(); i++) flux << ElementTable[i] 
                                                      << G4endl << G4endl;

   return flux;
}      
          
