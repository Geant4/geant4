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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#ifndef CellParameterisation_H
#define CellParameterisation_H 1

#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CellParameterisation : public G4VPVParameterisation
{ 
  public:
  
    CellParameterisation
    (G4Material* nucleus1, G4Material* cytoplasm1,
     G4Material* nucleus2, G4Material* cytoplasm2,
     G4Material* nucleus3, G4Material* cytoplasm3
     );

    virtual ~CellParameterisation();
   
    void ComputeTransformation (const G4int copyNo,G4VPhysicalVolume* physVol) const;
    
    void ComputeDimensions(G4Box&, 
				  const G4int, 
				  const G4VPhysicalVolume* ) const;

    void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Ellipsoid &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    G4int GetNoBoxes() const {return fPhantomTotalPixels;}

    G4Material* ComputeMaterial(const G4int copyNo,
                                      G4VPhysicalVolume* physVol,
                                const G4VTouchable*);

   // NEW
   
   G4int    GetPhantomTotalPixels()   const {return fPhantomTotalPixels;}  
   G4int    GetNucleusTotalPixels()   const {return fNucleusTotalPixels;}  
   G4int    GetCytoplasmTotalPixels() const {return fCytoplasmTotalPixels;}  
   G4double GetPixelSizeX()           const {return fDimCellBoxX;}  
   G4double GetPixelSizeY()           const {return fDimCellBoxY;}  
   G4double GetPixelSizeZ()           const {return fDimCellBoxZ;}  
   G4double GetCytoplasmMass()        const {return fCytoplasmMass;}  
   G4double GetNucleusMass()          const {return fNucleusMass;}  

   G4ThreeVector GetVoxelThreeVector(G4int i) const {return fMapCell[i];}
   G4double      GetMaterialVector(G4int i)   const {return fMaterial[i];}
   G4double      GetMassVector(G4int i)       const {return fMass[i];}
   G4int         GetTissueType(G4int i)       const {return fTissueType[i];}

   // SINGLETON

   static CellParameterisation * Instance()
   {
     return gInstance; 
   }
   
  private:

    G4Material * fNucleusMaterial1;
    G4Material * fCytoplasmMaterial1;
    G4Material * fNucleusMaterial2;
    G4Material * fCytoplasmMaterial2;
    G4Material * fNucleusMaterial3;
    G4Material * fCytoplasmMaterial3;
    
    G4VisAttributes * fNucleusAttributes1;
    G4VisAttributes * fCytoplasmAttributes1;
    G4VisAttributes * fNucleusAttributes2;
    G4VisAttributes * fCytoplasmAttributes2;
    G4VisAttributes * fNucleusAttributes3;
    G4VisAttributes * fCytoplasmAttributes3;
    
    G4ThreeVector * fMapCell;    // VOXEL COORDINATES
    G4double *      fMaterial;   // MATERIAL 
    G4double *      fMass;       // DENSITY REGION
    G4int *         fTissueType; // DENSITY REGION

    G4int fPhantomTotalPixels;
    G4int fNucleusTotalPixels;
    G4int fCytoplasmTotalPixels;
    
    G4double fDimCellBoxX;
    G4double fDimCellBoxY;
    G4double fDimCellBoxZ;
    G4double fNucleusMass;
    G4double fCytoplasmMass;
    
    static CellParameterisation* gInstance;
      
};

#endif


