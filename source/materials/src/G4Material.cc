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
// $Id: G4Material.cc 102843 2017-02-27 13:02:28Z gcosmo $
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
// 18-10-07, move definition of material index to InitialisePointers (V.Ivanchenko) 
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
  maxNbComponents        = fNumberOfComponents = fNumberOfElements = 1;
  fArrayLength           = maxNbComponents;
  theElementVector       = new G4ElementVector();

  const std::vector<G4String> elmnames = 
    G4NistManager::Instance()->GetNistElementNames();
  G4String enam, snam;
  G4int iz = G4lrint(z);
  if(iz < (G4int)elmnames.size()) { 
    snam = elmnames[iz]; 
    enam = snam; 
  } else { 
    enam = "ELM_" + name; 
    snam = name;
  }
  theElementVector->push_back(new G4Element(enam, snam, z, a));  

  fMassFractionVector    = new G4double[1];
  fMassFractionVector[0] = 1. ;
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
    
  maxNbComponents     = nComponents;
  fArrayLength        = maxNbComponents;
  fNumberOfComponents = fNumberOfElements = 0;
  theElementVector    = new G4ElementVector();
  theElementVector->reserve(maxNbComponents);  
    
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
  fChemicalFormula = fBaseMaterial->GetChemicalFormula();
  fMassOfMolecule  = fBaseMaterial->GetMassOfMolecule();

  fNumberOfElements = fBaseMaterial->GetNumberOfElements();     
  maxNbComponents = fNumberOfElements;
  fNumberOfComponents = fNumberOfElements;

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
    //delete fMaterialPropertiesTable;
    if(fMassFractionVector) { delete [] fMassFractionVector; } 
    if(fAtomsVector)        { delete [] fAtomsVector; } 
  }
  delete fIonisation; 
  if(VecNbOfAtomsPerVolume) { delete [] VecNbOfAtomsPerVolume; } 

  // Remove this material from theMaterialTable.
  //
  theMaterialTable[fIndexInTable] = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::InitializePointers()
{
  theElementVector         = nullptr;
  fMassFractionVector      = nullptr;
  fAtomsVector             = nullptr;
  fMaterialPropertiesTable = nullptr;
    
  VecNbOfAtomsPerVolume    = nullptr;
  fBaseMaterial            = nullptr;

  fChemicalFormula         = "";

  // initilized data members
  fDensity  = 0.0;
  fState    = kStateUndefined;
  fTemp     = 0.0;
  fPressure = 0.0;
  maxNbComponents     = 0;
  fArrayLength        = 0;
  fNumberOfComponents = 0;
  fNumberOfElements   = 0;
  TotNbOfAtomsPerVolume = 0.0;
  TotNbOfElectPerVolume = 0.0; 
  fRadlen = 0.0;
  fNuclInterLen = 0.0;
  fMassOfMolecule = 0.0;

  fIonisation = nullptr;
  fSandiaTable = nullptr;

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
  TotNbOfAtomsPerVolume = 0.;
  if (VecNbOfAtomsPerVolume) { delete [] VecNbOfAtomsPerVolume; }
  VecNbOfAtomsPerVolume = new G4double[fNumberOfElements];
  TotNbOfElectPerVolume = 0.;
  for (G4int i=0; i<fNumberOfElements; ++i) {
     Zi = (*theElementVector)[i]->GetZ();
     Ai = (*theElementVector)[i]->GetA();
     VecNbOfAtomsPerVolume[i] = Avogadro*fDensity*fMassFractionVector[i]/Ai;
     TotNbOfAtomsPerVolume += VecNbOfAtomsPerVolume[i];
     TotNbOfElectPerVolume += VecNbOfAtomsPerVolume[i]*Zi;
  }
        
  ComputeRadiationLength();
  ComputeNuclearInterLength();

  if (fIonisation) { delete fIonisation; }
  fIonisation  = new G4IonisParamMat(this);
  if (fSandiaTable) { delete fSandiaTable; }
  fSandiaTable = new G4SandiaTable(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::CopyPointersOfBaseMaterial()
{
  G4double factor = fDensity/fBaseMaterial->GetDensity();
  TotNbOfAtomsPerVolume = factor*fBaseMaterial->GetTotNbOfAtomsPerVolume();
  TotNbOfElectPerVolume = factor*fBaseMaterial->GetTotNbOfElectPerVolume();

  theElementVector = 
    const_cast<G4ElementVector*>(fBaseMaterial->GetElementVector());
  fMassFractionVector = 
    const_cast<G4double*>(fBaseMaterial->GetFractionVector());
  fAtomsVector = const_cast<G4int*>(fBaseMaterial->GetAtomsVector());

  const G4double* v = fBaseMaterial->GetVecNbOfAtomsPerVolume();
  if (VecNbOfAtomsPerVolume)  { delete [] VecNbOfAtomsPerVolume; }
  VecNbOfAtomsPerVolume = new G4double[fNumberOfElements];
  for (G4int i=0; i<fNumberOfElements; ++i) {
    VecNbOfAtomsPerVolume[i] = factor*v[i];
  }
  fRadlen = fBaseMaterial->GetRadlen()/factor;
  fNuclInterLen = fBaseMaterial->GetNuclearInterLength()/factor;

  if (fIonisation) { delete fIonisation; }
  fIonisation = new G4IonisParamMat(this);

  fSandiaTable = fBaseMaterial->GetSandiaTable();
  fIonisation->SetMeanExcitationEnergy(fBaseMaterial->GetIonisation()->GetMeanExcitationEnergy());

  fMaterialPropertiesTable = fBaseMaterial->GetMaterialPropertiesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// AddElement -- composition by atom count

void G4Material::AddElement(G4Element* element, G4int nAtoms)
{   
  // initialization
  if ( fNumberOfElements == 0 ) {
     fAtomsVector        = new G4int   [fArrayLength];
     fMassFractionVector = new G4double[fArrayLength];
  }

  // filling ...
  if ( fNumberOfElements < maxNbComponents ) {
     theElementVector->push_back(element);     
     fAtomsVector[fNumberOfElements] = nAtoms;
     fNumberOfComponents = ++fNumberOfElements;
  } else {
    G4cout << "G4Material::AddElement ERROR for " << fName << " nElement= " 
	   <<  fNumberOfElements << G4endl;
    G4Exception ("G4Material::AddElement()", "mat031", FatalException, 
           "Attempt to add more than the declared number of elements.");
  } 
  // filled.
  if ( fNumberOfElements == maxNbComponents ) {     
    // compute proportion by mass
    G4int i=0;
    G4double Amol = 0.;
    for (i=0; i<fNumberOfElements; ++i) {
      G4double w = fAtomsVector[i]*(*theElementVector)[i]->GetA(); 
      Amol += w;
      fMassFractionVector[i] = w;
    }
    for (i=0; i<fNumberOfElements; ++i) {
      fMassFractionVector[i] /= Amol;
    }

    fMassOfMolecule = Amol/Avogadro;
    ComputeDerivedQuantities();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// AddElement -- composition by fraction of mass

void G4Material::AddElement(G4Element* element, G4double fraction)
{
  if(fraction < 0.0 || fraction > 1.0) {
    G4cout << "G4Material::AddElement ERROR for " << fName << " and " 
	   << element->GetName() << "  mass fraction= " << fraction 
	   << " is wrong " << G4endl;
    G4Exception ("G4Material::AddElement()", "mat032", FatalException, 	   
                 "Attempt to add element with wrong mass fraction");
  }
  // initialization
  if (fNumberOfComponents == 0) {
    fMassFractionVector = new G4double[fArrayLength];
    fAtomsVector        = new G4int   [fArrayLength];
  }
  // filling ...
  if (fNumberOfComponents < maxNbComponents) {
    G4int el = 0;
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    while ((el<fNumberOfElements)&&(element!=(*theElementVector)[el])) { ++el; }
    if (el<fNumberOfElements) fMassFractionVector[el] += fraction;
    else {
      theElementVector->push_back(element); 
      fMassFractionVector[el] = fraction;
      ++fNumberOfElements;
    }
    ++fNumberOfComponents;  
  } else {
    G4cout << "G4Material::AddElement ERROR for " << fName << " nElement= " 
	   <<  fNumberOfElements << G4endl;
    G4Exception ("G4Material::AddElement()", "mat033", FatalException, 
           "Attempt to add more than the declared number of elements.");
  }    

  // filled.
  if (fNumberOfComponents == maxNbComponents) {

    G4int i=0;
    G4double Zmol(0.), Amol(0.);
    // check sum of weights -- OK?
    G4double wtSum(0.0);
    for (i=0; i<fNumberOfElements; ++i) {
      wtSum += fMassFractionVector[i];
      Zmol +=  fMassFractionVector[i]*(*theElementVector)[i]->GetZ();
      Amol +=  fMassFractionVector[i]*(*theElementVector)[i]->GetA();
    }
    if (std::fabs(1.-wtSum) > perThousand) {
      G4cerr << "WARNING !! for " << fName << " sum of fractional masses "
	     <<  wtSum << " is not 1 - results may be wrong" 
	     << G4endl;
    }
    for (i=0; i<fNumberOfElements; ++i) {
      fAtomsVector[i] = 
	G4lrint(fMassFractionVector[i]*Amol/(*theElementVector)[i]->GetA());
    }
     
    ComputeDerivedQuantities();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// AddMaterial -- composition by fraction of mass

void G4Material::AddMaterial(G4Material* material, G4double fraction)
{
  if(fraction < 0.0 || fraction > 1.0) {
    G4cout << "G4Material::AddMaterial ERROR for " << fName << " and " 
	   << material->GetName() << "  mass fraction= " << fraction 
	   << " is wrong ";
    G4Exception ("G4Material::AddMaterial()", "mat034", FatalException,
                 "Attempt to add material with wrong mass fraction");	   
  }
  // initialization
  if (fNumberOfComponents == 0) {
    fMassFractionVector = new G4double[fArrayLength];
    fAtomsVector        = new G4int   [fArrayLength];
  }

  G4int nelm = material->GetNumberOfElements();

  // arrays should be extended
  if(nelm > 1) {
    G4int nold    = fArrayLength;
    fArrayLength += nelm - 1;
    G4double* v1 = new G4double[fArrayLength];
    G4int* i1    = new G4int[fArrayLength];
    for(G4int i=0; i<nold; ++i) {
      v1[i] = fMassFractionVector[i];
      i1[i] = fAtomsVector[i];
    }
    delete [] fAtomsVector;
    delete [] fMassFractionVector;
    fMassFractionVector = v1;
    fAtomsVector = i1;
  }

  // filling ...
  if (fNumberOfComponents < maxNbComponents) {
    for (G4int elm=0; elm<nelm; ++elm)
      {
        G4Element* element = (*(material->GetElementVector()))[elm];
        G4int el = 0;
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
        while ((el<fNumberOfElements)&&(element!=(*theElementVector)[el])) el++;
        if (el < fNumberOfElements) fMassFractionVector[el] += fraction
                                          *(material->GetFractionVector())[elm];
        else {
	  theElementVector->push_back(element); 
          fMassFractionVector[el] = fraction
	                                  *(material->GetFractionVector())[elm];
          ++fNumberOfElements;
        }
      } 
    ++fNumberOfComponents;
    ///store massFraction of material component
    fMatComponents[material] = fraction;
      
  } else {
    G4cout << "G4Material::AddMaterial ERROR for " << fName << " nElement= " 
	   <<  fNumberOfElements << G4endl;
    G4Exception ("G4Material::AddMaterial()", "mat035", FatalException, 
           "Attempt to add more than the declared number of components.");
  }    

  // filled.
  if (fNumberOfComponents == maxNbComponents) {
    G4int i=0;
    G4double Zmol(0.), Amol(0.);
    // check sum of weights -- OK?
    G4double wtSum(0.0);
    for (i=0; i<fNumberOfElements; ++i) {
      wtSum += fMassFractionVector[i];
      Zmol +=  fMassFractionVector[i]*(*theElementVector)[i]->GetZ();
      Amol +=  fMassFractionVector[i]*(*theElementVector)[i]->GetA();
    }
    if (std::fabs(1.-wtSum) > perThousand) {
      G4cout << "G4Material::AddMaterial WARNING !! for " << fName 
	     << " sum of fractional masses "
	     <<  wtSum << " is not 1 - results may be wrong" 
	     << G4endl;
    }
    for (i=0; i<fNumberOfElements; ++i) {
      fAtomsVector[i] = 
	G4lrint(fMassFractionVector[i]*Amol/(*theElementVector)[i]->GetA());
    }
     
    ComputeDerivedQuantities();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeRadiationLength()
{
  G4double radinv = 0.0 ;
  for (G4int i=0;i<fNumberOfElements;++i) {
     radinv += VecNbOfAtomsPerVolume[i]*((*theElementVector)[i]->GetfRadTsai());
  }
  fRadlen = (radinv <= 0.0 ? DBL_MAX : 1./radinv);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4Material::ComputeNuclearInterLength()
{
  static const G4double lambda0  = 35*CLHEP::g/CLHEP::cm2;
  static const G4double twothird = 2.0/3.0;
  G4double NILinv = 0.0;
  for (G4int i=0; i<fNumberOfElements; ++i) {
    G4int Z = G4lrint( (*theElementVector)[i]->GetZ());
    G4double A = (*theElementVector)[i]->GetN();
    if(1 == Z) {
      NILinv += VecNbOfAtomsPerVolume[i]*A;
    } else {
      NILinv += VecNbOfAtomsPerVolume[i]*G4Exp(twothird*G4Log(A));
    } 
  }
  NILinv *= amu/lambda0; 
  fNuclInterLen = (NILinv <= 0.0 ? DBL_MAX : 1./NILinv);
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

G4Material* 
G4Material::GetMaterial(const G4String& materialName, G4bool warning)
{  
  // search the material by its name 
  for (size_t J=0 ; J<theMaterialTable.size() ; ++J)
    {
      if (theMaterialTable[J]->GetName() == materialName)
	{ return theMaterialTable[J]; }
    }
   
  // the material does not exist in the table
  if (warning) {
    G4cout << "G4Material::GetMaterial() WARNING: The material: "
	   << materialName 
	   << " does not exist in the table. Return NULL pointer."
	   << G4endl;
  }	 
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4Material::GetZ() const
{ 
  if (fNumberOfElements > 1) {
     G4cout << "G4Material ERROR in GetZ. The material: " << fName 
	    << " is a mixture.";
     G4Exception ("G4Material::GetZ()", "mat036", FatalException, 
                  "the Atomic number is not well defined." );
  } 
  return (*theElementVector)[0]->GetZ();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4Material::GetA() const
{ 
  if (fNumberOfElements > 1) { 
    G4cout << "G4Material ERROR in GetA. The material: " << fName 
	   << " is a mixture.";
    G4Exception ("G4Material::GetA()", "mat037", FatalException,  
		 "the Atomic mass is not well defined." );
  } 
  return  (*theElementVector)[0]->GetA();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4ExtendedMaterial.hh"

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
    << G4BestUnit(material->GetIonisation()->GetMeanExcitationEnergy(),
		  "Energy")
    << "  temperature: " << std::setw(6) << std::setprecision(2)  
    << (material->fTemp)/kelvin << " K"
    << "  pressure: "    << std::setw(6) << std::setprecision(2)   
    << (material->fPressure)/atmosphere << " atm" << "\n";
  
  for (G4int i=0; i<material->fNumberOfElements; i++) {
    flux 
      << "\n   ---> " << (*(material->theElementVector))[i] 
      << "\n          ElmMassFraction: " 
      << std::setw(6)<< std::setprecision(2) 
      << (material->fMassFractionVector[i])/perCent << " %" 
      << "  ElmAbundance "     << std::setw(6)<< std::setprecision(2) 
      << 100*(material->VecNbOfAtomsPerVolume[i])
      /(material->TotNbOfAtomsPerVolume)
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
     
std::ostream& operator<<(std::ostream& flux, G4MaterialTable MaterialTable)
{
  //Dump info for all known materials
  flux << "\n***** Table : Nb of materials = " << MaterialTable.size() 
       << " *****\n" << G4endl;
        
  for (size_t i=0; i<MaterialTable.size(); ++i) { 
    flux << MaterialTable[i] << G4endl << G4endl; 
  }

  return flux;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4Material::IsExtended() const
{ return false; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
