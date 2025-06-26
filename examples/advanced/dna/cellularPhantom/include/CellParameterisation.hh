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
// -----------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet (a), S. Incerti (a), N. H. Tran (a), L. Morelli (a,b)
//
// a) University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
// b) Politecnico di Milano, Italy
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// -----------------------------------------------------------------------------

#ifndef CellParameterisation_H
#define CellParameterisation_H 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CellParameterisation : public G4VPVParameterisation
{
  public:

    explicit CellParameterisation
      (G4String fileName,
       G4Material *RedMat, G4Material *GreenMat, G4Material *BlueMat,
       G4double shiftX, G4double shiftY, G4double shiftZ);

    ~CellParameterisation() override;

    void ComputeTransformation
      (const G4int copyNo, G4VPhysicalVolume *physVol) const override;

    G4Material *ComputeMaterial (const G4int copyNo,
                                       G4VPhysicalVolume *physVol,
                                 const G4VTouchable *) override;

    inline auto GetPhantomTotalPixels() const { return fPhantomTotalPixels; }

    inline auto GetRedTotalPixels() const { return fRedTotalPixels; }
    inline auto GetGreenTotalPixels() const { return fGreenTotalPixels; }
    inline auto GetBlueTotalPixels() const { return fBlueTotalPixels; }

    inline auto GetPixelSizeX() const { return fDimCellBoxX; }
    inline auto GetPixelSizeY() const { return fDimCellBoxY; }
    inline auto GetPixelSizeZ() const { return fDimCellBoxZ; }

    inline auto GetRedMass() const { return fRedMass; }
    inline auto GetGreenMass() const { return fGreenMass; }
    inline auto GetBlueMass() const { return fBlueMass; }

    //inline auto GetVoxelThreeVector(G4int i) const { return fMapCell[i]; }
    //inline auto GetVoxelThreeVectorPixel(G4int i) const { return fMapCellPxl[i]; }
    inline auto GetVoxelThreeVectorOriginal(G4int i) const { return fMapCellOriginal[i]; }

    inline auto GetMaterial(G4int i) const { return fMaterial[i]; }

    // Singleton
    static CellParameterisation *Instance()
    {
      return gInstance;
    }

  private:

    void Initialize(const G4String&);

    static CellParameterisation *gInstance;

    G4double fDimCellBoxX = 0;
    G4double fDimCellBoxY = 0;
    G4double fDimCellBoxZ = 0;

    G4double fSizeRealX = 0;
    G4double fSizeRealY = 0;
    G4double fSizeRealZ = 0;

    G4Material *fRedMaterial = nullptr;
    G4Material *fGreenMaterial = nullptr;
    G4Material *fBlueMaterial = nullptr;

    G4double fShiftX = 0.;
    G4double fShiftY = 0.;
    G4double fShiftZ = 0.;

    G4VisAttributes *fRedAttributes = nullptr;
    G4VisAttributes *fGreenAttributes = nullptr;
    G4VisAttributes *fBlueAttributes = nullptr;

    G4ThreeVector *fMapCell = nullptr; // VOXEL COORDINATES
    G4ThreeVector *fMapCellPxl = nullptr;// VOXEL COORDINATES IN PIXEL, NO SHIFT
    G4ThreeVector *fMapCellOriginal = nullptr; // VOXEL COORDINATES (original space)

    G4int    *fMaterial = nullptr; // MATERIAL

    G4int fPhantomTotalPixels = 0;
    G4int fRedTotalPixels = 0;
    G4int fGreenTotalPixels = 0;
    G4int fBlueTotalPixels = 0;

    G4double fRedMass = 0.;
    G4double fGreenMass = 0.;
    G4double fBlueMass = 0.;

    char fRealUnit;

    G4double fOffsetX = 0.;
    G4double fOffsetY = 0.;
    G4double fOffsetZ = 0.;
};

#endif
