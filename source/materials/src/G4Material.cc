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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// 26-06-96, Code uses operators (+=, *=, ++, -> etc.) correctly, P. Urban
// 10-07-96, new data members added by L.Urban
// 12-12-96, new data members added by L.Urban
// 20-01-97, aesthetic rearrangement. RadLength calculation modified.
//           Data members Zeff and Aeff REMOVED (i.e. passed to the Elements).
//           (local definition of Zeff in DensityEffect and FluctModel...)
//           Vacuum defined as a G4State. Mixture flag removed, M.Maire. 
// 29-01-97, State=Vacuum automatically set density=0 in the contructors.
//           Subsequent protections have been put in the calculation of 
//           MeanExcEnergy, ShellCorrectionVector, DensityEffect, M.Maire.
// 11-02-97, ComputeDensityEffect() rearranged, M.Maire.
// 20-03-97, corrected initialization of pointers, M.Maire.
// 28-05-98, the kState=kVacuum has been removed.
//           automatic check for a minimal density, M.Maire 
// 12-06-98, new method AddMaterial() allowing mixture of materials, M.Maire  
// 09-07-98, ionisation parameters removed from the class, M.Maire
// 05-10-98, change names: NumDensity -> NbOfAtomsPerVolume
// 18-11-98, new interface to SandiaTable
// 19-01-99  enlarge tolerance on test of coherence of gas conditions
// 19-07-99, Constructors with chemicalFormula added by V.Ivanchenko
// 16-01-01, Nuclear interaction length, M.Maire
// 12-03-01, G4bool fImplicitElement;
//           copy constructor and assignement operator revised (mma)
// 03-05-01, flux.precision(prec) at begin/end of operator<<
// 17-07-01, migration to STL. M. Verderi.
// 14-09-01, Suppression of the data member fIndexInTable
// 26-02-02, fIndexInTable renewed
// 16-04-02, G4Exception put in constructor with chemical formula
// 06-05-02, remove the check of the ideal gas state equation
// 06-08-02, remove constructors with chemical formula (mma)
// 22-01-04, proper STL handling of theElementVector (Hisaya)
// 30-03-05, warning in GetMaterial(materialName) 
// 09-03-06, minor change of printout (V.Ivanchenko) 
// 10-01-07, compute fAtomVector in the case of mass fraction (V.Ivanchenko) 
// 27-07-07, improve destructor (V.Ivanchenko) 
// 18-10-07, moved definition of mat index to InitialisePointers (V.Ivanchenko) 
// 13-08-08, do not use fixed size arrays (V.Ivanchenko)
// 26-10-11, new scheme for G4Exception  (mma)
// 13-04-12, map<G4Material*,G4double> fMatComponents, filled in AddMaterial()
// 21-04-12, fMassOfMolecule, computed for AtomsCount (mma)
// 
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iomanip>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ExtendedMaterial.hh"
#include "G4AtomicShells.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

G4MaterialTable G4Material::theMaterialTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to create a material from scratch

G4Material::G4Material(const G4String& name, G4double z,
                       G4double a, G4double density, 
                       G4State state, G4double temp, G4double pressure)
  : fName(name)		       
{
  InitializePointers();
    
  if (density < universe_mean_density)
    { 
      G4cout << " G4Material WARNING:"
	     << " define a material with density=0 is not allowed. \n"
	     << " The material " << name << " will be constructed with the"
	     << " default minimal density: " << universe_mean_density/(g/cm3) 
	     << "g/cm3" << G4endl;
      density = universe_mean_density;
    } 

  fDensity  = density;
  fState    = state;
  fTemp     = temp;
  fPressure = pressure;

  // Initialize theElementVector allocating one
  // element corresponding to this material
  fNbComponents = fNumberOfElements = 1;
  theElementVector = new G4ElementVector();

  // take element from DB
  G4NistManager* nist = G4NistManager::Instance();
  G4int iz = G4lrint(z);
  auto elm = nist->FindOrBuildElement(iz);
  if(elm == nullptr)
  {
    elm = new G4Element("ELM_" + name, name, z, a);
  }
  theElementVector->push_back(elm);  

  fMassFractionVector    = new G4double[1];
  fMassFractionVector[0] = 1.;
  fMassOfMolecule        = a/CLHEP::Avogadro;
  
  if (fState == kStateUndefined)
    {
      if (fDensity > kGasThreshold) { fState = kStateSolid; }
      else                          { fState = kStateGas; }
    }

  ComputeDerivedQuantities();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to create a material from a List of constituents
// (elements and/or materials)  added with AddElement or AddMaterial

G4Material::G4Material(const G4String& name, G4double density,
                       G4int nComponents,
                       G4State state, G4double temp, G4double pressure)
  : fName(name)		       
{
  InitializePointers();
    
  if (density < universe_mean_density)
    {
      G4cout << "--- Warning from G4Material::G4Material()"
	     << " define a material with density=0 is not allowed. \n"
	     << " The material " << name << " will be constructed with the"
	     << " default minimal density: " << universe_mean_density/(g/cm3) 
	     << "g/cm3" << G4endl;
      density = universe_mean_density;
    }
        
  fDensity  = density;
  fState    = state;
  fTemp     = temp;
  fPressure = pressure;
    
  fNbComponents = nComponents;
  fMassFraction = true;
    
  if (fState == kStateUndefined) 
    {
      if (fDensity > kGasThreshold) { fState = kStateSolid; }
      else                          { fState = kStateGas; }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor to create a material from base material

G4Material::G4Material(const G4String& name, G4double density,
                       const G4Material* bmat,
                       G4State state, G4double temp, G4double pressure)
  : fName(name)		       
{
  InitializePointers();
    
  if (density < universe_mean_density)
    {
      G4cout << "--- Warning from G4Material::G4Material()"
	     << " define a material with density=0 is not allowed. \n"
	     << " The material " << name << " will be constructed with the"
	     << " default minimal density: " << universe_mean_density/(g/cm3) 
	     << "g/cm3" << G4endl;
      density = universe_mean_density;
    }

  fDensity  = density;
  fState    = state;
  fTemp     = temp;
  fPressure = pressure;

  fBaseMaterial = bmat;
  auto ptr = bmat;
  if(nullptr != ptr) {
    while(1) {
      ptr = ptr->GetBaseMaterial();
      if(nullptr == ptr) { break; }
      else { fBaseMaterial = ptr; }
    }
  }

  fChemicalFormula = fBaseMaterial->GetChemicalFormula();
  fMassOfMolecule  = fBaseMaterial->GetMassOfMolecule();

  fNumberOfElements = (G4int)fBaseMaterial->GetNumberOfElements();     
  fNbComponents = fNumberOfElements;

  CopyPointersOfBaseMaterial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Material::G4Material(__void__&)
  : fName("")
{
  InitializePointers();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material::~G4Material()
{
  //  G4cout << "### Destruction of material " << fName << " started" <<G4endl;
  if(fBaseMaterial == nullptr) {
    delete theElementVector;  
    delete fSandiaTable; 
    delete [] fMassFractionVector;
    delete [] fAtomsVector;
  }
  delete fIonisation; 
  delete [] fVecNbOfAtomsPerVolume; 

  // Remove this material from theMaterialTable.
  //
  theMaterialTable[fIndexInTable] = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::InitializePointers()
{
  fBaseMaterial            = nullptr;
  fMaterialPropertiesTable = nullptr;
  theElementVector         = nullptr;
  fAtomsVector             = nullptr;
  fMassFractionVector      = nullptr;
  fVecNbOfAtomsPerVolume   = nullptr;

  fIonisation  = nullptr;
  fSandiaTable = nullptr;

  fDensity = fFreeElecDensity = fTemp = fPressure = 0.0;
  fTotNbOfAtomsPerVolume = 0.0;
  fTotNbOfElectPerVolume = 0.0; 
  fRadlen = fNuclInterLen = fMassOfMolecule = 0.0;
    
  fState = kStateUndefined;
  fNumberOfElements = fNbComponents = fIdxComponent = 0;
  fMassFraction = true;
  fChemicalFormula = "";

  // Store in the static Table of Materials
  fIndexInTable = theMaterialTable.size();
  for(size_t i=0; i<fIndexInTable; ++i) {
    if(theMaterialTable[i]->GetName() == fName) {
      G4cout << "G4Material WARNING: duplicate name of material "
	     << fName << G4endl;
      break;
    }
  }
  theMaterialTable.push_back(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeDerivedQuantities()
{
  // Header routine to compute various properties of material.
  // 
  // Number of atoms per volume (per element), total nb of electrons per volume
  G4double Zi, Ai;
  fTotNbOfAtomsPerVolume = 0.;
  delete[] fVecNbOfAtomsPerVolume;
  fVecNbOfAtomsPerVolume = new G4double[fNumberOfElements];
  fTotNbOfElectPerVolume = 0.;
  fFreeElecDensity = 0.;
  const G4double elecTh = 15.*CLHEP::eV; // threshold for conductivity e-
  for (G4int i=0; i<fNumberOfElements; ++i) {
    Zi = (*theElementVector)[i]->GetZ();
    Ai = (*theElementVector)[i]->GetA();
    fVecNbOfAtomsPerVolume[i] = Avogadro*fDensity*fMassFractionVector[i]/Ai;
    fTotNbOfAtomsPerVolume += fVecNbOfAtomsPerVolume[i];
    fTotNbOfElectPerVolume += fVecNbOfAtomsPerVolume[i]*Zi;
    if(fState != kStateGas) {
      fFreeElecDensity += fVecNbOfAtomsPerVolume[i]*
	G4AtomicShells::GetNumberOfFreeElectrons(Zi, elecTh);
    }
  }
        
  ComputeRadiationLength();
  ComputeNuclearInterLength();

  if(fIonisation == nullptr)
  {
    fIonisation = new G4IonisParamMat(this);
  }
  if(fSandiaTable == nullptr)
  {
    fSandiaTable = new G4SandiaTable(this);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::CopyPointersOfBaseMaterial()
{
  G4double factor = fDensity/fBaseMaterial->GetDensity();
  fTotNbOfAtomsPerVolume = factor*fBaseMaterial->GetTotNbOfAtomsPerVolume();
  fTotNbOfElectPerVolume = factor*fBaseMaterial->GetTotNbOfElectPerVolume();
  fFreeElecDensity = factor*fBaseMaterial->GetFreeElectronDensity();

  if(fState == kStateUndefined) { fState = fBaseMaterial->GetState(); }

  theElementVector =
    const_cast<G4ElementVector*>(fBaseMaterial->GetElementVector());
  fMassFractionVector = 
    const_cast<G4double*>(fBaseMaterial->GetFractionVector());
  fAtomsVector = const_cast<G4int*>(fBaseMaterial->GetAtomsVector());

  const G4double* v = fBaseMaterial->GetVecNbOfAtomsPerVolume();
  delete[] fVecNbOfAtomsPerVolume;
  fVecNbOfAtomsPerVolume = new G4double[fNumberOfElements];
  for (G4int i=0; i<fNumberOfElements; ++i) {
    fVecNbOfAtomsPerVolume[i] = factor*v[i];
  }
  fRadlen = fBaseMaterial->GetRadlen()/factor;
  fNuclInterLen = fBaseMaterial->GetNuclearInterLength()/factor;

  if(fIonisation == nullptr)
  {
    fIonisation = new G4IonisParamMat(this);
  }
  fIonisation->SetMeanExcitationEnergy(fBaseMaterial->GetIonisation()->GetMeanExcitationEnergy());
  if(fBaseMaterial->GetIonisation()->GetDensityEffectCalculator() != nullptr)
  {
    ComputeDensityEffectOnFly(true);
  }

  fSandiaTable = fBaseMaterial->GetSandiaTable();
  fMaterialPropertiesTable = fBaseMaterial->GetMaterialPropertiesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4Material::AddElementByNumberOfAtoms(const G4Element* elm, G4int nAtoms)
{   
  // perform checks consistency
  if(0 == fIdxComponent) {
    fMassFraction = false;
    fAtoms = new std::vector<G4int>;
    fElm = new std::vector<const G4Element*>;
  }
  if(fIdxComponent >= fNbComponents) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element "
       << elm->GetName() << " with Natoms=" << nAtoms
       << " wrong attempt to add more than the declared number of elements "
       << fIdxComponent << " >= " << fNbComponents;
    G4Exception ("G4Material::AddElementByNumberOfAtoms()", "mat031",
                 FatalException, ed, "");
  }
  if(fMassFraction) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element "
       << elm->GetName() << " with Natoms=" << nAtoms
       << " problem: cannot add by number of atoms after "
       << "addition of elements by mass fraction";
    G4Exception ("G4Material::AddElementByNumberOfAtoms()", "mat031",
                 FatalException, ed, "");
  }
  if(0 >= nAtoms) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element "
       << elm->GetName() << " with Natoms=" << nAtoms
       << " problem: number of atoms should be above zero";
    G4Exception ("G4Material::AddElementByNumberOfAtoms()", "mat031",
                 FatalException, ed, "");
  }

  // filling
  G4bool isAdded = false;
  if(!fElm->empty()) {
    for (G4int i=0; i<fNumberOfElements; ++i) {
      if ( elm == (*fElm)[i] ) {
	(*fAtoms)[i] += nAtoms;
	isAdded = true;
	break;
      }
    }
  }
  if(!isAdded) {
    fElm->push_back(elm);     
    fAtoms->push_back(nAtoms);
    ++fNumberOfElements;
  } 
  ++fIdxComponent;

  // is filled - complete composition of atoms
  if (fIdxComponent == fNbComponents) {     
    theElementVector = new G4ElementVector();
    theElementVector->reserve(fNumberOfElements);
    fAtomsVector = new G4int[fNumberOfElements];
    fMassFractionVector = new G4double[fNumberOfElements];
    
    G4double Amol = 0.;
    for (G4int i=0; i<fNumberOfElements; ++i) {
      theElementVector->push_back((*fElm)[i]);
      fAtomsVector[i] = (*fAtoms)[i];
      G4double w = fAtomsVector[i]*(*fElm)[i]->GetA(); 
      Amol += w;
      fMassFractionVector[i] = w;
    }
    for (G4int i=0; i<fNumberOfElements; ++i) {
      fMassFractionVector[i] /= Amol;
    }
    delete fAtoms;
    delete fElm;
    fMassOfMolecule = Amol/CLHEP::Avogadro;
    ComputeDerivedQuantities();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4Material::AddElementByMassFraction(const G4Element* elm, G4double fraction)
{
  // perform checks consistency
  if(fraction < 0.0 || fraction > 1.0) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element "
       << elm->GetName() << " massFraction= " << fraction
       << " is wrong ";
    G4Exception ("G4Material::AddElementByMassFraction()", "mat031",
		 FatalException, ed, "");
  }
  if(!fMassFraction) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element "
       << elm->GetName() << ", massFraction= " << fraction
       << ", fIdxComponent=" << fIdxComponent
       << " problem: cannot add by mass fraction after "
       << "addition of elements by number of atoms";
    G4Exception ("G4Material::AddElementByMassFraction()", "mat031",
		 FatalException, ed, "");
  }
  if(fIdxComponent >= fNbComponents) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added element " 
       << elm->GetName() << ", massFraction= " << fraction
       << ", fIdxComponent=" << fIdxComponent
       << "; attempt to add more than the declared number of components "
       << fIdxComponent << " >= " << fNbComponents;
    G4Exception ("G4Material::AddElementByMassFraction()", "mat031",
		 FatalException, ed, "");
  }
  if(0 == fIdxComponent) {
    fElmFrac = new std::vector<G4double>;
    fElm = new std::vector<const G4Element*>;
  }

  // filling 
  G4bool isAdded = false;
  if(!fElm->empty()) {
    for (G4int i=0; i<fNumberOfElements; ++i) {
      if ( elm == (*fElm)[i] ) {
	(*fElmFrac)[i] += fraction;
	isAdded = true;
        break;
      }
    }
  }
  if(!isAdded) {
    fElm->push_back(elm);
    fElmFrac->push_back(fraction);
    ++fNumberOfElements;
  }
  ++fIdxComponent;

  // is filled
  if(fIdxComponent == fNbComponents) { FillVectors(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// composition by fraction of mass
void G4Material::AddMaterial(G4Material* material, G4double fraction)
{
  if(fraction < 0.0 || fraction > 1.0) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added material " 
       << material->GetName() << ", massFraction= " << fraction
       << " is wrong ";
    G4Exception ("G4Material::AddMaterial()", "mat031", FatalException,
                 ed, "");	   
  }
  if(!fMassFraction) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added material " 
       << material->GetName() << ", massFraction= " << fraction
       << ", fIdxComponent=" << fIdxComponent
       << " problem: cannot add by mass fraction after "
       << "addition of elements by number of atoms";
    G4Exception ("G4Material::AddMaterial()", "mat031", FatalException,
		 ed, "");
  }
  if(fIdxComponent >= fNbComponents) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " and added material " 
       << material->GetName() << ", massFraction= " << fraction
       << "; attempt to add more than the declared number of components "
       << fIdxComponent << " >= " << fNbComponents;
    G4Exception ("G4Material::AddMaterial()", "mat031", FatalException, 
		 ed, "");
  }
  if(0 == fIdxComponent) {
    fElmFrac = new std::vector<G4double>;
    fElm = new std::vector<const G4Element*>;
  }

  // filling 
  G4int nelm = (G4int)material->GetNumberOfElements();
  for(G4int j=0; j<nelm; ++j) {
    auto elm = material->GetElement(j);
    auto frac = material->GetFractionVector();
    G4bool isAdded = false;
    if(!fElm->empty()) {
      for (G4int i=0; i<fNumberOfElements; ++i) {
	if ( elm == (*fElm)[i] ) {
	  (*fElmFrac)[i] += fraction*frac[j];
	  isAdded = true;
	  break;
	}
      }
    }
    if(!isAdded) {
      fElm->push_back(elm);     
      fElmFrac->push_back(fraction*frac[j]);
      ++fNumberOfElements;
    }
  }

  fMatComponents[material] = fraction;
  ++fIdxComponent;

  // is filled
  if(fIdxComponent == fNbComponents) { FillVectors(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::FillVectors()
{
  // there are material components
  theElementVector = new G4ElementVector();
  theElementVector->reserve(fNumberOfElements);
  fAtomsVector = new G4int[fNumberOfElements];
  fMassFractionVector = new G4double[fNumberOfElements];
    
  G4double wtSum(0.0);
  for (G4int i=0; i<fNumberOfElements; ++i) {
    theElementVector->push_back((*fElm)[i]);
    fMassFractionVector[i] = (*fElmFrac)[i];
    wtSum += fMassFractionVector[i];
  }
  delete fElmFrac;
  delete fElm;

  // check sum of weights -- OK?
  if (std::abs(1.-wtSum) > perThousand) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " sum of fractional masses "
	   <<  wtSum << " is not 1 - results may be wrong";
    G4Exception ("G4Material::FillVectors()", "mat031", JustWarning, 
		 ed, "");
  }
  G4double coeff = (wtSum > 0.0) ? 1./wtSum : 1.0;
  G4double Amol(0.);
  for (G4int i=0; i<fNumberOfElements; ++i) {
    fMassFractionVector[i] *= coeff;
    Amol += fMassFractionVector[i]*(*theElementVector)[i]->GetA();
  }
  for (G4int i=0; i<fNumberOfElements; ++i) {
    fAtomsVector[i] = 
      G4lrint(fMassFractionVector[i]*Amol/(*theElementVector)[i]->GetA());
  }
  ComputeDerivedQuantities();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeRadiationLength()
{
  G4double radinv = 0.0 ;
  for (G4int i=0; i<fNumberOfElements; ++i) {
    radinv += fVecNbOfAtomsPerVolume[i]*((*theElementVector)[i]->GetfRadTsai());
  }
  fRadlen = (radinv <= 0.0 ? DBL_MAX : 1./radinv);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeNuclearInterLength()
{
  const G4double lambda0  = 35*CLHEP::g/CLHEP::cm2;
  const G4double twothird = 2.0/3.0;
  G4double NILinv = 0.0;
  for (G4int i=0; i<fNumberOfElements; ++i) {
    G4int Z = (*theElementVector)[i]->GetZasInt();
    G4double A = (*theElementVector)[i]->GetN();
    if(1 == Z) {
      NILinv += fVecNbOfAtomsPerVolume[i]*A;
    } else {
      NILinv += fVecNbOfAtomsPerVolume[i]*G4Exp(twothird*G4Log(A));
    } 
  }
  NILinv *= amu/lambda0; 
  fNuclInterLen = (NILinv <= 0.0 ? DBL_MAX : 1./NILinv);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::SetChemicalFormula(const G4String& chF)  
{
  if(!IsLocked()) { fChemicalFormula = chF; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::SetFreeElectronDensity(G4double val)
{
  if(val >= 0. && !IsLocked()) { fFreeElecDensity = val; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeDensityEffectOnFly(G4bool val)
{
  if(!IsLocked()) {
    if (nullptr == fIonisation) {
      fIonisation  = new G4IonisParamMat(this);
    }
    fIonisation->ComputeDensityEffectOnFly(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MaterialTable* G4Material::GetMaterialTable()
{
  return &theMaterialTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t G4Material::GetNumberOfMaterials()
{
  return theMaterialTable.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4Material::GetMaterial(const G4String& materialName, G4bool warn)
{  
  // search the material by its name
  for(auto& j : theMaterialTable)
  {
    if(j->GetName() == materialName)
    {
      return j;
    }
  }

  // the material does not exist in the table
  if (warn) {
    G4cout << "G4Material::GetMaterial() WARNING: The material: "
	   << materialName 
	   << " does not exist in the table. Return NULL pointer."
	   << G4endl;
  }	 
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4Material::GetMaterial(G4double z, G4double a, G4double dens)
{  
  // search the material by its name
  for(auto mat : theMaterialTable)
  {
    if(1 == mat->GetNumberOfElements() && z == mat->GetZ() &&
       a == mat->GetA() && dens == mat->GetDensity())
    {
      return mat;
    }
  }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4Material::GetMaterial(size_t nComp, G4double dens)
{  
  // search the material by its name
  for(auto mat : theMaterialTable)
  {
    if(nComp == mat->GetNumberOfElements() && dens == mat->GetDensity())
    {
      return mat;
    }
  }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4Material::GetZ() const
{ 
  if (fNumberOfElements > 1) {
    G4ExceptionDescription ed;
    ed << "For material " << fName << " ERROR in GetZ() - Nelm=" 
       << fNumberOfElements << " > 1, which is not allowed"; 
    G4Exception ("G4Material::GetZ()", "mat036", FatalException, 
		 ed, "");
  } 
  return (*theElementVector)[0]->GetZ();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4Material::GetA() const
{ 
  if (fNumberOfElements > 1) { 
    G4ExceptionDescription ed;
    ed << "For material " << fName << " ERROR in GetA() - Nelm=" 
       << fNumberOfElements << " > 1, which is not allowed"; 
    G4Exception ("G4Material::GetA()", "mat036", FatalException, 
		 ed, "");
  } 
  return (*theElementVector)[0]->GetA();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream& operator<<(std::ostream& flux, const G4Material* material)
{
  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  G4long prec = flux.precision(3);
  
  flux
    << " Material: "         << std::setw(8) <<  material->fName
    << " " << material->fChemicalFormula << " "
    << "  density: "         << std::setw(6) << std::setprecision(3)  
    << G4BestUnit(material->fDensity,"Volumic Mass") 
    << "  RadL: "            << std::setw(7)  << std::setprecision(3)  
    << G4BestUnit(material->fRadlen,"Length")
    << "  Nucl.Int.Length: " << std::setw(7)  << std::setprecision(3)  
    << G4BestUnit(material->fNuclInterLen,"Length") 
    << "\n" << std::setw(30)   
    << "  Imean: "           << std::setw(7)  << std::setprecision(3)  
    << G4BestUnit(material->GetIonisation()->GetMeanExcitationEnergy(),"Energy")
    << "  temperature: " << std::setw(6) << std::setprecision(2)  
    << (material->fTemp)/CLHEP::kelvin << " K"
    << "  pressure: "    << std::setw(6) << std::setprecision(2)   
    << (material->fPressure)/CLHEP::atmosphere << " atm" << "\n";
  
  for (G4int i=0; i<material->fNumberOfElements; i++) {
    flux 
      << "\n   ---> " << (*(material->theElementVector))[i] 
      << "\n          ElmMassFraction: " 
      << std::setw(6)<< std::setprecision(2) 
      << (material->fMassFractionVector[i])/perCent << " %" 
      << "  ElmAbundance "     << std::setw(6)<< std::setprecision(2) 
      << 100*(material->fVecNbOfAtomsPerVolume[i])
      /(material->fTotNbOfAtomsPerVolume)
      << " % \n";
  }
  flux.precision(prec);    
  flux.setf(mode,std::ios::floatfield);

  if(material->IsExtended())
  { static_cast<const G4ExtendedMaterial*>(material)->Print(flux); }

  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::ostream& operator<<(std::ostream& flux, const G4Material& material)
{
  flux << &material;        
  return flux;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
     
std::ostream& operator<<(std::ostream& flux, const G4MaterialTable& MaterialTable)
{
  //Dump info for all known materials
  flux << "\n***** Table : Nb of materials = " << MaterialTable.size() 
       << " *****\n" << G4endl;

  for(auto i : MaterialTable)
  {
    flux << i << G4endl << G4endl;
  }

  return flux;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4Material::IsExtended() const
{ 
  return false; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::SetMaterialPropertiesTable(G4MaterialPropertiesTable* anMPT)
{
  if(fMaterialPropertiesTable != anMPT && !IsLocked()) {
    delete fMaterialPropertiesTable;
    fMaterialPropertiesTable = anMPT;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4Material::IsLocked()
{
  auto state = G4StateManager::GetStateManager()->GetCurrentState();
  return !(state == G4State_PreInit ||
	   state == G4State_Init || 
	   state == G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
