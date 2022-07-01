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

//---------------------------------------------------------------------------
//
// ClassName:   G4CrystalExtension
//
// Description: Contains crystal properties
//
// Class description:
//
// Extension of G4Material for the management of a crystal
// structure. It has to be attached to a G4ExtendedMaterial
// in order to instantiate a G4LogicalCrystalVolume.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 21-04-16, created by E.Bagli

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4CrystalExtension_HH
#define G4CrystalExtension_HH 1

#include "G4VMaterialExtension.hh"
#include "G4CrystalAtomBase.hh"
#include "G4AtomicBond.hh"
#include "G4CrystalUnitCell.hh"
#include "G4NistManager.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4CrystalExtension : public G4VMaterialExtension
{
public:  // with description
    //
    // Constructor to create a material
    //
    G4CrystalExtension(G4Material* ,const G4String& name = "crystal");

    ~G4CrystalExtension() override;

    void Print() const override { ; };

   private:
    G4Material* fMaterial;
    
public:
    G4Material* GetMaterial() {return fMaterial;};
    void SetMaterial(G4Material* aMat) {fMaterial = aMat;};
    
    
    //
    // Crystal cell description, i.e. space group
    // and cell parameters
    //

private:
    G4CrystalUnitCell* theUnitCell;
public:
    inline void SetUnitCell(G4CrystalUnitCell* aUC) {theUnitCell=aUC;}
    inline G4CrystalUnitCell* GetUnitCell()
        const {return theUnitCell;}
    
    //
    // Elasticity and reduced elasticity tensors
    //
    
public:
    typedef G4double Elasticity[3][3][3][3];
    typedef G4double ReducedElasticity[6][6];
    
protected:
    Elasticity fElasticity;	    	    // Full 4D elasticity tensor
    ReducedElasticity fElReduced;		    // Reduced 2D elasticity tensor
    
public:
    const Elasticity& GetElasticity() const { return fElasticity; }
    const ReducedElasticity& GetElReduced() const { return fElReduced; }
    
public:
    G4double GetCijkl(G4int i, G4int j, G4int k, G4int l) const {
        return fElasticity[i][j][k][l];
    }
    
    // Reduced elasticity tensor: C11-C66 interface for clarity
    void SetElReduced(const ReducedElasticity& mat);
    
    void SetCpq(G4int p, G4int q, G4double value);
    G4double GetCpq(G4int p, G4int q) const { return fElReduced[p-1][q-1]; }

    //
    // Map of atom positions for each element
    // The coordinate system is the unit cell
    //
private:
    std::map<const G4Element*,G4CrystalAtomBase*> theCrystalAtomBaseMap;

public:
    G4CrystalAtomBase* GetAtomBase(const G4Element* anElement);
    
    void AddAtomBase(const G4Element* anElement,G4CrystalAtomBase* aBase){
        theCrystalAtomBaseMap.insert(std::pair<const G4Element*,G4CrystalAtomBase*>(anElement,aBase));
    }


    G4CrystalAtomBase* GetAtomBase(G4int anElIdx){
        return GetAtomBase(fMaterial->GetElement(anElIdx));
    }

    void AddAtomBase(G4int anElIdx,G4CrystalAtomBase* aLattice){
        AddAtomBase(fMaterial->GetElement(anElIdx),aLattice);
    }

    //
    // Get the position of all the atoms in the unit cell
    // for a single element or all the elements
    //
    G4bool GetAtomPos(const G4Element* anEl, std::vector<G4ThreeVector>& vecout);
    G4bool GetAtomPos(std::vector<G4ThreeVector>& vecout);
    
    G4bool GetAtomPos(G4int anElIdx, std::vector<G4ThreeVector>& vecout){
        GetAtomPos(fMaterial->GetElement(anElIdx),vecout);
        return true;
    }

    
    //
    // Structure factor calculations
    // Eq. 46, Chapter 2 , Introduction to solid state physics, C. Kittel
    //
    G4complex ComputeStructureFactor(G4double kScatteringVector,
                                     G4int h,
                                     G4int k,
                                     G4int l);
    G4complex ComputeStructureFactorGeometrical(G4int h,
                                                G4int k,
                                                G4int l);
    
    //
    // Bond between atoms. Each bond is mapped with two Elements
    // and the number of the atoms in the corresponding G4CrystalBaseAtomPos
    //
    
private:
    std::vector<G4AtomicBond*> theAtomicBondVector;

public:
    void AddAtomicBond(G4AtomicBond* aBond) {theAtomicBondVector.push_back(aBond);};
    G4AtomicBond* GetAtomicBond(G4int idx) {return theAtomicBondVector[idx];};
    std::vector<G4AtomicBond*> GetAtomicBondVector() {return theAtomicBondVector;};
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
