// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Material.hh,v 1.4 1999-07-28 11:48:18 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//      ------------------- class G4Material ---------------------
//
//                     Torre Wenaus, November 1995
//
// Materials defined via the G4Material class are used to define the
// composition of Geant volumes.
// a Material is always made of Elements. It can be defined directly
// from scratch (defined by a single element), specifying :
//                                             its name,
//                                             density,
//                                             state informations,
//                                             and Z,A of the underlying Element.
// or in terms of a collection of constituent Elements with specified weights
// (composition specified either by fractional mass or atom counts).
//
// Quantities, with physical meaning or not, which are constant in a given 
// material are computed and stored here as Derived data members.
// The class contains as a private static member the table of defined
// materials (an ordered vector of materials).

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4MATERIAL_HH
#define G4MATERIAL_HH

#include "G4ios.hh"
#include <rw/tpvector.h>
#include <rw/tpordvec.h>
#include "globals.hh"
#include "G4Element.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4IonisParamMat.hh"
#include "G4SandiaTable.hh"

typedef RWTPtrVector<G4Element> G4ElementVector;

class G4Material;              //forward declaration
typedef RWTPtrOrderedVector<G4Material> G4MaterialTable;

enum G4State { kStateUndefined, kStateSolid, kStateLiquid, kStateGas };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4Material
{
public:

    //
    // Constructor to create a material from scratch.
    //
    G4Material(const G4String& name,				//its name
                     G4double  z, 				//atomic number
                     G4double  a,				//mass of mole
                     G4double  density, 			//density
                     G4State   state    = kStateUndefined,	//solid,liqid,gas
                     G4double  temp     = STP_Temperature,	//temperature
                     G4double  pressure = STP_Pressure);	//pressure

    //
    // Constructor to create a material from a combination of elements
    // and/or materials subsequently added via AddElement and/or AddMaterial
    //
    G4Material(const G4String& name,				//its name
                     G4double  density, 			//density
                     G4int     nComponents,			//nb of components 
                     G4State   state    = kStateUndefined,	//solid,liquid,gas
                     G4double  temp     = STP_Temperature,	//temperature
                     G4double  pressure = STP_Pressure);	//pressure

    //
    // Constructor to create a material with chemical formula from scratch.
    //
    G4Material(const G4String& name,				//its name
	       const G4String& chFormula,                       //chemical formula
                     G4double  z, 				//atomic number
                     G4double  a,				//mass of mole
                     G4double  density, 			//density
                     G4State   state    = kStateUndefined,	//solid,liqid,gas
                     G4double  temp     = STP_Temperature,	//temperature
                     G4double  pressure = STP_Pressure);	//pressure

    //
    // Constructor to create a material with chemical formula from a 
    // combination of elements and/or materials subsequently added via 
    // AddElement and/or AddMaterial
    //
    G4Material(const G4String& name,				//its name
	       const G4String& chFormula,                       //chemical formula
                     G4double  density, 			//density
                     G4int     nComponents,			//nb of components 
                     G4State   state    = kStateUndefined,	//solid,liquid,gas
                     G4double  temp     = STP_Temperature,	//temperature
                     G4double  pressure = STP_Pressure);	//pressure

    //
    // Add an element, giving number of atoms
    //
    void AddElement(G4Element* element,				//the element
                    G4int      nAtoms);				//nb of atoms in a molecule

    //
    // Add an element or material, giving fraction of mass
    //
    void AddElement (G4Element* element ,			//the element
                     G4double   fraction);			//fraction of mass
                     
    void AddMaterial(G4Material* material,			//the material
                     G4double   fraction);			//fraction of mass
                     
                     
   ~G4Material();
                        
    //
    // retrieval methods
    // 
    G4String GetName()            const {return fName;};
    G4String GetChemicalFormula() const {return fChemicalFormula;};
    G4double GetDensity()         const {return fDensity;};

    G4State  GetState()       const {return fState;};
    G4double GetTemperature() const {return fTemp;};
    G4double GetPressure()    const {return fPressure;};
    
    const
    G4ElementVector* GetElementVector()  const {return theElementVector;};
    size_t GetNumberOfElements()         const {return fNumberOfElements;};  
    const  G4double* GetFractionVector() const {return fMassFractionVector;};
    const  G4int*    GetAtomsVector()    const {return fAtomsVector;};

    const 
    G4Element* GetElement(G4int iel) const {return (*theElementVector)[iel];};

    void SetMaterialPropertiesTable(G4MaterialPropertiesTable* anMPT)
                                       {fMaterialPropertiesTable = anMPT;};
    				       
    G4MaterialPropertiesTable* GetMaterialPropertiesTable() const
                                       {return fMaterialPropertiesTable;};

    static
    const G4MaterialTable* GetMaterialTable() {return &theMaterialTable;};
    
    static
    size_t GetNumberOfMaterials()   {return theMaterialTable.length();};
    size_t GetIndex() const         {return fIndexInTable;};
    
    static  G4Material* GetMaterial(G4String name);

    const
    G4double* GetVecNbOfAtomsPerVolume() const {return VecNbOfAtomsPerVolume;};
    G4double  GetTotNbOfAtomsPerVolume() const {return TotNbOfAtomsPerVolume;};
    G4double  GetTotNbOfElectPerVolume() const {return TotNbOfElectPerVolume;};

    //obsolete names (5-10-98) see the 2 functions above
    const
    G4double* GetAtomicNumDensityVector() const {return VecNbOfAtomsPerVolume;};
    G4double  GetElectronDensity()        const {return TotNbOfElectPerVolume;};
    
    G4double         GetRadlen()          const {return fRadlen;};
    G4IonisParamMat* GetIonisation()      const {return fIonisation;};
    G4SandiaTable*   GetSandiaTable()     const {return fSandiaTable;};
    
    G4double GetZ() const;
    G4double GetA() const;
    
    friend
    ostream& operator<<(ostream&, G4Material*);
    
    friend
    ostream& operator<<(ostream&, G4Material&);
    
    friend
    ostream& operator<<(ostream&, G4MaterialTable);
        
    G4int operator==(const G4Material&) const;
    G4int operator!=(const G4Material&) const;

private:

    G4Material(const G4Material&);
    const G4Material& operator=(const G4Material&);

    void InitializePointers();
     
    // Header routine for all derived quantities
    void ComputeDerivedQuantities();

    // Compute Radiation length
    void ComputeRadiationLength();

private:

//
// Basic data members ( To define a material)
//

    G4String         fName;                 // Material name
    G4String         fChemicalFormula;      // Material chemical formula
    G4double         fDensity;              // Material density
   
    G4State          fState;                // Material state (defaults to undefined,  
                                            //   determined internally based on density)
    G4double         fTemp;                 // Temperature (defaults to STP)
    G4double         fPressure;             // Pressure    (defaults to STP)

    G4int            maxNbComponents;       // total number of components in the material 
    size_t           fNumberOfComponents;   // Number of components declared so far

    size_t           fNumberOfElements;     // Number of Elements in the material
    G4ElementVector* theElementVector;      // vector of constituent Elements
    G4double*        fMassFractionVector;   // composition by fractional mass
    G4int*           fAtomsVector;          // composition by atom count

    G4MaterialPropertiesTable* fMaterialPropertiesTable;

    static
    G4MaterialTable theMaterialTable;       // the material table
    size_t          fIndexInTable;          // Index of material in the material table

//
// Derived data members (computed from the basic data members)
//
    // some general atomic properties
    
    G4double* VecNbOfAtomsPerVolume;      // vector of nb of atoms per volume
    G4double  TotNbOfAtomsPerVolume;      // total nb of atoms per volume 
    G4double  TotNbOfElectPerVolume;      // total nb of electrons per volume 
    G4double  fRadlen;                    // Radiation length  
    
    G4IonisParamMat* fIonisation;          // ionisation parameters
    G4SandiaTable*   fSandiaTable;         // Sandia table         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4Material* G4Material::GetMaterial(G4String materialName)
{  
  // search the material by its name 
  for (G4int J=0 ; J<theMaterialTable.length() ; J++)
   {
    if(theMaterialTable[J]->GetName() == materialName)
      return theMaterialTable[J];
   }
   
  G4cerr << "  Warning from GetMaterial(name). The material: " << materialName
         << "  does not exist in the MaterialTable.  Return NULL pointer \n";
  return NULL;          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4double G4Material::GetZ() const
{ 
  if (fNumberOfElements > 1) {
     G4cerr << "WARNING in GetZ. The material: " << fName << " is a mixture." << endl;
     G4Exception ( " the Atomic number is not well defined." );
  } 
  return (*theElementVector)(0)->GetZ();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....


inline
G4double G4Material::GetA() const
{ 
  if (fNumberOfElements > 1) { 
     G4cerr << "WARNING in GetA. The material: " << fName << " is a mixture." << endl;
     G4Exception ( " the Atomic mass is not well defined." );
  } 
  return  (*theElementVector)(0)->GetA();      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#endif
