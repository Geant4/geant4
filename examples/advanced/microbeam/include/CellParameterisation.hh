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
  
    explicit CellParameterisation
    (G4Material* nucleus1, G4Material* cytoplasm1,
     G4Material* nucleus2, G4Material* cytoplasm2,
     G4Material* nucleus3, G4Material* cytoplasm3
     );

    ~CellParameterisation() override;
   
    void ComputeTransformation (const G4int copyNo,G4VPhysicalVolume* physVol) const override;
    
    void ComputeDimensions(G4Box&, 
				  const G4int, 
				  const G4VPhysicalVolume* ) const override;

    void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override {}
	
    void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override {}
	
    void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Ellipsoid &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const override{}

    G4int GetNoBoxes() {return fPhantomTotalPixels;}

    G4Material* ComputeMaterial(const G4int copyNo,
                                      G4VPhysicalVolume* physVol,
                                const G4VTouchable*) override;

   // NEW
   
   G4int   GetPhantomTotalPixels()   {return fPhantomTotalPixels;}  
   G4int   GetNucleusTotalPixels()   {return fNucleusTotalPixels;}  
   G4int   GetCytoplasmTotalPixels() {return fCytoplasmTotalPixels;}  
   G4double GetPixelSizeX() {return fDimCellBoxX;}  
   G4double GetPixelSizeY() {return fDimCellBoxY;}  
   G4double GetPixelSizeZ() {return fDimCellBoxZ;}  
   G4double GetCytoplasmMass() {return fCytoplasmMass;}  
   G4double GetNucleusMass()   {return fNucleusMass;}  

   G4ThreeVector GetVoxelThreeVector(G4int i) {return fMapCell[i];}
   G4double GetMaterialVector(G4int i) {return fMaterial[i];}
   G4double GetMassVector(G4int i) {return fMass[i];}
   G4int GetTissueType(G4int i) {return fTissueType[i];}

//SINGLETON

   static CellParameterisation * Instance()
   {
     return gInstance; 
   }
//
   
  private:
    static CellParameterisation* gInstance;
    
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
    G4double fDimCellBoxX;
    G4double fDimCellBoxY;
    G4double fDimCellBoxZ;
    G4double fNucleusMass;
    G4double fCytoplasmMass;
    
    G4int *         fTissueType; // DENSITY REGION
    G4int fPhantomTotalPixels;
    G4int fNucleusTotalPixels;
    G4int fCytoplasmTotalPixels;      
};
#endif


