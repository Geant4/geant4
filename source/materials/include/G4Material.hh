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
//---------------------------------------------------------------------------
//
// ClassName:   G4Material
//
// Description: Contains material properties
//
// Class description:
//
// Is used to define the material composition of Geant4 volumes.
// A G4Material is always made of G4Elements. It should has the name, 
// the list of G4Elements, material density, material state, temperature, 
// pressure. Other parameters are optional and may be set by the user code 
// or computed at initialisation. 
// 
// There is several ways to construct G4Material:
//   - from single element;
//   - from a list of components (elements or other materials);
//   - from internal Geant4 database of materials
//
// A collection of constituent Elements/Materials should be defined 
// with specified weights by fractional mass or atom counts (only for Elements).
//
// Quantities, with physical meaning or not, which are constant in a given 
// material are computed and stored here as Derived data members.
//
// The class contains as a private static member the Table of defined
// materials (an ordered vector of materials).
//
// It is strongly not recommended to delete materials in user code.
// All materials will be deleted automatically at the end of Geant4 session.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 10-07-96, new data members added by L.Urban
// 12-12-96, new data members added by L.Urban
// 20-01-97, aesthetic rearrangement. RadLength calculation modified
//           Data members Zeff and Aeff REMOVED (i.e. passed to the Elements).
//           (local definition of Zeff in DensityEffect and FluctModel...)
//           Vacuum defined as a G4State. Mixture flag removed, M.Maire  
// 29-01-97, State=Vacuum automatically set density=0 in the contructors.
//           Subsequent protections have been put in the calculation of 
//           MeanExcEnergy, ShellCorrectionVector, DensityEffect, M.Maire
// 20-03-97, corrected initialization of pointers, M.Maire
// 10-06-97, new data member added by V.Grichine (fSandiaPhotoAbsCof)
// 27-06-97, new function GetElement(int), M.Maire
// 24-02-98, fFractionVector become fMassFractionVector
// 28-05-98, kState=kVacuum removed: 
//           The vacuum is an ordinary gas vith very low density, M.Maire
// 12-06-98, new method AddMaterial() allowing mixture of materials, M.Maire
// 09-07-98, Ionisation parameters removed from the class, M.Maire
// 04-08-98, new method GetMaterial(materialName), M.Maire
// 05-10-98, change name: NumDensity -> NbOfAtomsPerVolume
// 18-11-98, SandiaTable interface modified.
// 19-07-99, new data member (chemicalFormula) added by V.Ivanchenko
// 12-03-01, G4bool fImplicitElement (mma)
// 30-03-01, suppression of the warning message in GetMaterial
// 17-07-01, migration to STL. M. Verderi.
// 14-09-01, Suppression of the data member fIndexInTable
// 31-10-01, new function SetChemicalFormula() (mma)
// 26-02-02, fIndexInTable renewed
// 06-08-02, remove constructors with ChemicalFormula (mma)
// 15-11-05, GetMaterial(materialName, G4bool warning=true)
// 13-04-12, std::map<G4Material*,G4double> fMatComponents (mma)
// 21-04-12, fMassOfMolecule (mma)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MATERIAL_HH
#define G4MATERIAL_HH 1

#include <vector>
#include <map>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ios.hh"
#include "G4Element.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4IonisParamMat.hh"
#include "G4SandiaTable.hh"
#include "G4ElementVector.hh"
#include "G4MaterialTable.hh"

enum G4State { kStateUndefined = 0, kStateSolid, kStateLiquid, kStateGas };

static const G4double NTP_Temperature = 293.15*CLHEP::kelvin;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Material
{
public:  // with description
  //
  // Constructor to create a material from single element
  //
  G4Material(const G4String& name,				//its name
                   G4double  z, 				//atomic number
                   G4double  a,					//mass of mole
                   G4double  density, 				//density
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

  //
  // Constructor to create a material from a combination of elements
  // and/or materials subsequently added via AddElement and/or AddMaterial
  //
  G4Material(const G4String& name,				//its name
                   G4double  density, 				//density
                   G4int     nComponents,			//nbOfComponents
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

  //
  // Constructor to create a material from the base material
  //
  G4Material(const G4String& name,				//its name
                   G4double  density, 				//density
             const G4Material* baseMaterial,			//base material
                   G4State   state    = kStateUndefined,	//solid,gas
                   G4double  temp     = NTP_Temperature,	//temperature
                   G4double  pressure = CLHEP::STP_Pressure);	//pressure

  //
  // Add an element, giving number of atoms
  //
  void AddElementByNumberOfAtoms(const G4Element* elm, G4int nAtoms); 
  inline 
  void AddElement(G4Element* elm, G4int nAtoms)
  { AddElementByNumberOfAtoms(elm, nAtoms); }

  //
  // Add an element or material, giving fraction of mass
  //
  void AddElementByMassFraction(const G4Element* elm, G4double fraction);  
  inline void AddElement (G4Element* elm, G4double frac)
  { AddElementByMassFraction(elm, frac); }

  void AddMaterial(G4Material* material, G4double fraction);
                     
  virtual ~G4Material();  
  //
  // retrieval methods
  // 
  inline const G4String& GetName()            const {return fName;}
  inline const G4String& GetChemicalFormula() const {return fChemicalFormula;}
  inline G4double GetFreeElectronDensity()    const {return fFreeElecDensity;}
  inline G4double GetDensity()     const {return fDensity;}
  inline G4State  GetState()       const {return fState;}
  inline G4double GetTemperature() const {return fTemp;}
  inline G4double GetPressure()    const {return fPressure;}
    
  //number of elements constituing this material:    
  inline size_t GetNumberOfElements()  const {return fNumberOfElements;}
    
  //vector of pointers to elements constituing this material:          
  inline const
  G4ElementVector* GetElementVector()  const {return theElementVector;}
  
  //vector of fractional mass of each element:
  inline const  
  G4double* GetFractionVector() const {return fMassFractionVector;}
    
  //vector of atom count of each element:
  inline const  
  G4int*    GetAtomsVector()    const {return fAtomsVector;}

  //return a pointer to an element, given its index in the material:
  inline const 
  G4Element* GetElement(G4int iel) const {return (*theElementVector)[iel];}
  
  //vector of nb of atoms per volume of each element in this material:
  inline const
  G4double* GetVecNbOfAtomsPerVolume() const {return fVecNbOfAtomsPerVolume;}
  //total number of atoms per volume:
  inline
  G4double  GetTotNbOfAtomsPerVolume() const {return fTotNbOfAtomsPerVolume;}
  //total number of electrons per volume:
  inline
  G4double  GetTotNbOfElectPerVolume() const {return fTotNbOfElectPerVolume;}

  //obsolete names (5-10-98) see the 2 functions above
  inline const
  G4double* GetAtomicNumDensityVector() const {return fVecNbOfAtomsPerVolume;}
  inline G4double  GetElectronDensity() const {return fTotNbOfElectPerVolume;}
    
  // Radiation length:     
  inline G4double  GetRadlen()            const {return fRadlen;}
    
  // Nuclear interaction length     
  inline G4double GetNuclearInterLength() const {return fNuclInterLen;}
        
  // ionisation parameters:
  inline G4IonisParamMat* GetIonisation() const {return fIonisation;}

  // Sandia table:
  inline G4SandiaTable* GetSandiaTable()  const {return fSandiaTable; }
  
  // Base material:
  inline 
  const G4Material* GetBaseMaterial()     const {return fBaseMaterial;}
  
  // material components:
  inline
  const std::map<G4Material*,G4double>& GetMatComponents() const 
                                                {return fMatComponents;}
					       
  // for chemical compound
  inline G4double GetMassOfMolecule() const     {return fMassOfMolecule;}

  void SetChemicalFormula(const G4String& chF);

  void SetFreeElectronDensity(G4double val);

  void ComputeDensityEffectOnFly(G4bool);
      
  // meaningful only for single material:
  G4double GetZ() const;
  G4double GetA() const;

  //the MaterialPropertiesTable (if any) attached to this material:
  void SetMaterialPropertiesTable(G4MaterialPropertiesTable* anMPT);
  				       
  inline G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
  {return fMaterialPropertiesTable;}

  //the index of this material in the Table:    
  inline size_t GetIndex() const {return fIndexInTable;}

  // the static Table of Materials:
  //
  static G4MaterialTable* GetMaterialTable();
      
  static size_t GetNumberOfMaterials();
      
  //return  pointer to a material, given its name:    
  static G4Material* GetMaterial(const G4String& name, G4bool warning=true);

  //return  pointer to a simple material, given its propeties:
  static G4Material* GetMaterial(G4double z, G4double a, G4double dens);

  //return  pointer to a composit material, given its propeties:
  static G4Material* GetMaterial(size_t nComp, G4double dens);
  
  //
  //printing methods
  //
  friend std::ostream& operator<<(std::ostream&, const G4Material*);    
  friend std::ostream& operator<<(std::ostream&, const G4Material&);    
  friend std::ostream& operator<<(std::ostream&, const G4MaterialTable&);
    
  G4Material(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  inline void SetName (const G4String& name) {fName=name;}

  virtual G4bool IsExtended() const;

  // operators
  G4bool operator==(const G4Material&) const = delete;
  G4bool operator!=(const G4Material&) const = delete;
  G4Material(const G4Material&) = delete;
  const G4Material& operator=(const G4Material&) = delete;

private:

  void InitializePointers();
   
  // Header routine for all derived quantities
  void ComputeDerivedQuantities();

  // Compute Radiation length
  void ComputeRadiationLength();
  
  // Compute Nuclear interaction length
  void ComputeNuclearInterLength();

  // Copy pointers of base material
  void CopyPointersOfBaseMaterial();

  void FillVectors();

  G4bool IsLocked();

  static
  G4MaterialTable theMaterialTable;  // the material table

  const G4Material* fBaseMaterial;   // Pointer to the base material
  G4MaterialPropertiesTable* fMaterialPropertiesTable;

  //
  // General atomic properties defined in constructor or
  // computed from the basic data members
  //
  
  G4ElementVector* theElementVector;// vector of constituent G4Elements
  G4int* fAtomsVector;             // composition by atom count
  G4double* fMassFractionVector;   // composition by fractional mass   
  G4double* fVecNbOfAtomsPerVolume;// number of atoms per volume
  
  G4IonisParamMat* fIonisation;    // ionisation parameters
  G4SandiaTable* fSandiaTable;     // Sandia table         

  G4double fDensity;               // Material density
  G4double fFreeElecDensity;       // Free electron density
  G4double fTemp;                  // Temperature (defaults: STP)
  G4double fPressure;              // Pressure    (defaults: STP)

  G4double fTotNbOfAtomsPerVolume; // Total nb of atoms per volume 
  G4double fTotNbOfElectPerVolume; // Total nb of electrons per volume 
  G4double fRadlen;                // Radiation length
  G4double fNuclInterLen;          // Nuclear interaction length  
  G4double fMassOfMolecule;        // Correct for materials built by atoms count

  G4State fState;                  // Material state
  size_t fIndexInTable;            // Index in the material table 
  G4int fNumberOfElements;         // Number of G4Elements in the material

  // Class members used only at initialisation
  G4int fNbComponents;             // Number of components 
  G4int fIdxComponent;             // Index of a new component
  G4bool fMassFraction;            // Flag of the method to add components

  // For composites built 
  std::vector<G4int>* fAtoms = nullptr;
  std::vector<G4double>* fElmFrac = nullptr;
  std::vector<const G4Element*>* fElm = nullptr;

  // For composites built via AddMaterial()
  std::map<G4Material*, G4double> fMatComponents; 

  G4String fName;                  // Material name
  G4String fChemicalFormula;       // Material chemical formula
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
