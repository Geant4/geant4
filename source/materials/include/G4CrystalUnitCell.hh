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
//
//---------------------------------------------------------------
//
// G4CrystalUnitCell
//
// Class Description:
//


#ifndef G4CrystalUnitCell_H
#define G4CrystalUnitCell_H 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "G4CrystalBravaisLattices.h"
#include "G4CrystalLatticeSystems.h"
//#include "sginfo.h"

class G4CrystalUnitCell
{
public:
    G4CrystalUnitCell(G4double sizeA,
                      G4double sizeB,
                      G4double sizeC,
                      G4double alpha,
                      G4double beta,
                      G4double gamma,
                      G4int spacegroup);

    virtual ~G4CrystalUnitCell();
    
private:
    G4int theSpaceGroup; //
public:
    inline G4int GetSpaceGroup() const {return theSpaceGroup;};
    inline void SetSpaceGroup(G4int aInt)  {theSpaceGroup=aInt;};

private:
    theLatticeSystemType GetLatticeSystem(G4int aGroup);
    theBravaisLatticeType GetBravaisLattice(G4int aGroup);

public:
    theLatticeSystemType GetLatticeSystem(){
        return GetLatticeSystem(theSpaceGroup);
    }
    theBravaisLatticeType GetBravaisLattice(){
        return GetBravaisLattice(theSpaceGroup);
    }

private:
    //T_SgInfo SgInfo;
    /*!< struct from SgInfo library needed for further calculations
     see SgInfo documentation on http://cci.lbl.gov/sginfo/ */
private:
    G4double cosa,cosb,cosg;
    G4double sina,sinb,sing;
    G4double cosar,cosbr,cosgr;
    
    //
    // Size and angles of the crystalline unit cell
    //
protected:
    G4ThreeVector nullVec;
    
    G4ThreeVector theSize; // cell sizes
    G4ThreeVector theAngle; // cell angles
    G4ThreeVector theUnitBasis[3];	// Basis unit vectors in direct orientation
    G4ThreeVector theBasis[3];	// Basis vectors in direct orientation

public:
    const G4ThreeVector& GetBasis(G4int idx) const;
    const G4ThreeVector& GetUnitBasis(G4int idx) const;
    inline G4ThreeVector GetSize() const {return theSize;}
    inline G4ThreeVector GetAngle() const {return theAngle;}
   
    G4ThreeVector GetUnitBasisTrigonal(); // return theUnitBase[2] vector
    //
    // Reciprocal size and angles of the crystalline unit cell
    //
protected:
    G4ThreeVector theRecSize; // reciprocal cell sizes
    G4ThreeVector theRecAngle; // reciprocal cell angles
    G4ThreeVector theRecUnitBasis[3];	// Basis unit vectors in reciprocal orientation
    G4ThreeVector theRecBasis[3];	// Basis vectors in reciprocal orientation
    
public:
    const G4ThreeVector& GetRecBasis(G4int idx) const;
    const G4ThreeVector& GetRecUnitBasis(G4int idx) const;
    inline G4ThreeVector GetRecSize() const {return theRecSize;}
    inline G4ThreeVector GetRecAngle() const {return theRecAngle;}

    //
    // Methods to populate atom position in the lattice from the basis
    // and the unit basis
    //
public:
    G4bool FillAtomicUnitPos(G4ThreeVector& pos, std::vector<G4ThreeVector>& vecout);
    G4bool FillAtomicPos(G4ThreeVector& pos, std::vector<G4ThreeVector>& vecout);

    //
    // Methods to populate elasticity and reduced elasticity tensors
    //
public:
    G4bool FillElReduced(G4double Cij[6][6]);
private:
    G4bool FillAmorphous(G4double Cij[6][6]) const;
    G4bool FillCubic(G4double Cij[6][6]) const;
    G4bool FillTetragonal(G4double Cij[6][6]) const;
    G4bool FillOrthorhombic(G4double Cij[6][6]) const;
    G4bool FillRhombohedral(G4double Cij[6][6]) const;
    G4bool FillMonoclinic(G4double Cij[6][6]) const;
    G4bool FillTriclinic(G4double Cij[6][6]) const;
    G4bool FillHexagonal(G4double Cij[6][6]) const;
    
    G4bool ReflectElReduced(G4double Cij[6][6]) const;

    //
    // The volumes of the cell
    //
public:
    G4double ComputeCellVolume(); //compute and store the volume

    inline G4double GetVolume() const {return theVolume;} //get the stored volume
    inline G4double GetRecVolume() const {return theRecVolume;} //get the stored volume
    
private:
    G4double theVolume; // the cell volume
    G4double theRecVolume; // the cell volume
    
    //
    // Squared Reciprocal and direct interplanar spacing
    //
public:
    G4double GetIntSp2(G4int h,
                       G4int k,
                       G4int l); // squared interplanar spacing
    
    G4double GetRecIntSp2(G4int h,
                          G4int k,
                          G4int l); // squared reciprocal interplanar spacing
    
    G4double GetIntCosAng(G4int h1,
                          G4int k1,
                          G4int l1,
                          G4int h2,
                          G4int k2,
                          G4int l2); // cosine of the angle between two planes
    
};

#endif

