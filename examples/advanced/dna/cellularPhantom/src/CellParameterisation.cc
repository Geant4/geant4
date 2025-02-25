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
// --------------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet, S. Incerti, N. H. Tran, L. Morelli
//
// University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// --------------------------------------------------------------------------------

#include "CellParameterisation.hh"

#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation *CellParameterisation::gInstance = nullptr;

CellParameterisation::CellParameterisation
(G4String fileName,
 G4Material *RedMat, G4Material *GreenMat, G4Material *BlueMat,
 G4double shiftX, G4double shiftY, G4double shiftZ
)
:fRedMaterial(RedMat), fGreenMaterial(GreenMat), fBlueMaterial(BlueMat),
 fShiftX(shiftX), fShiftY(shiftY), fShiftZ(shiftZ)
{
  Initialize(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::Initialize(const G4String &fileName)
{
  G4int ncols, l, mat;
  G4int pixelX, pixelY, pixelZ;
  G4double x, y, z, den1, den2, den3;

  ncols = 0;
  l = 0;

  // Read phantom

  FILE *fMap;
  fMap = fopen(fileName, "r");

  fRedMass = 0;
  fGreenMass = 0;
  fBlueMass = 0;

  ncols = fscanf(fMap, "%d  %d  %d  %d", &fPhantomTotalPixels, &fRedTotalPixels, &fGreenTotalPixels,
                                         &fBlueTotalPixels);
  ncols = fscanf(fMap, "%lf  %lf  %lf  %s", &fSizeRealX, &fSizeRealY, &fSizeRealZ, &fRealUnit);
  ncols = fscanf(fMap, "%lf  %lf  %lf  %s", &fDimCellBoxX, &fDimCellBoxY, &fDimCellBoxZ, &fRealUnit);

  fMapCell = new G4ThreeVector[fPhantomTotalPixels]; //geant4 coordinates space
  fMapCellPxl = new G4ThreeVector[fPhantomTotalPixels]; //voxel space
  fMapCellOriginal = new G4ThreeVector[fPhantomTotalPixels]; //original coordinates space
  fMaterial = new G4int[fPhantomTotalPixels];

  fDimCellBoxX = fDimCellBoxX * um;
  fDimCellBoxY = fDimCellBoxY * um;
  fDimCellBoxZ = fDimCellBoxZ * um;

  den1 = fRedMaterial->GetDensity();
  den2 = fGreenMaterial->GetDensity();
  den3 = fBlueMaterial->GetDensity();

  fOffsetX = -fSizeRealX / 2 *um;
  fOffsetY = -fSizeRealY / 2 *um;
  fOffsetZ = -fSizeRealZ / 2 *um;

  G4cout << G4endl;
  G4cout << " #########################################################################" << G4endl;
  G4cout << "                               Phantom placement and density              " << G4endl;
  G4cout << " #########################################################################" << G4endl;
  G4cout << G4endl;
  G4cout << " ==========> Phantom origin - X (um) = " << (fOffsetX + fShiftX)/um << G4endl;
  G4cout << " ==========> Phantom origin - Y (um) = " << (fOffsetY + fShiftY)/um << G4endl;
  G4cout << " ==========> Phantom origin - Z (um) = " << (fOffsetZ + fShiftZ)/um << G4endl;
  G4cout << G4endl;
  G4cout << " ==========> Red density (g/cm3)   = " << den1/(g/cm3) << G4endl;
  G4cout << " ==========> Green density (g/cm3) = " << den2/(g/cm3) << G4endl;
  G4cout << " ==========> Blue density (g/cm3)  = " << den3/(g/cm3) << G4endl;
  G4cout << G4endl;
  G4cout << " #########################################################################" << G4endl;
  G4cout << G4endl;

  while (1)
  {
    ncols = fscanf(fMap, "%lf %lf %lf %d", &x, &y, &z, &mat);
    if (ncols < 0) break;

    G4ThreeVector v( x*um + fOffsetX + fShiftX, // phantom shift
                   -(y*um + fOffsetY + fShiftY),
                     z*um + fOffsetZ + fShiftZ );

    // Pixel coordinates
    pixelX = (x*um)/fDimCellBoxX;
    pixelY = (y*um)/fDimCellBoxY;
    pixelZ = (z*um)/fDimCellBoxZ;

    G4ThreeVector w(pixelX, pixelY, pixelZ);

    G4ThreeVector v_original(x*um, y*um, z*um);

    fMapCell[l] = v;
    fMapCellPxl[l] = w;
    fMapCellOriginal[l] = v_original;

    fMaterial[l] = mat;

    if (mat == 1){
      fRedMass += den1 * fDimCellBoxX * fDimCellBoxY * fDimCellBoxZ;
    }
    else if (mat == 2){
      fGreenMass += den2 * fDimCellBoxX * fDimCellBoxY * fDimCellBoxZ;
    }
    else if (mat == 3){
      fBlueMass += den3 * fDimCellBoxX * fDimCellBoxY * fDimCellBoxZ;
    }
    l++;
  }

  fclose(fMap);

  fRedAttributes = new G4VisAttributes;
  fRedAttributes->SetColour(G4Colour(1, 0, 0));
  fRedAttributes->SetForceSolid(false);

  fGreenAttributes = new G4VisAttributes;
  fGreenAttributes->SetColour(G4Colour(0, 1, 0));
  fGreenAttributes->SetForceSolid(false);

  fBlueAttributes = new G4VisAttributes;
  fBlueAttributes->SetColour(G4Colour(0, 0, 1));
  fBlueAttributes->SetForceSolid(false);

  gInstance = this;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::~CellParameterisation()
{
  delete[] fMapCell;
  delete[] fMapCellPxl;
  delete[] fMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume *physVol) const
{
  if(fMapCell == nullptr)
  {
    G4ExceptionDescription ex;
    ex<< "fMapCell == nullptr ";
    G4Exception("CellParameterisation::ComputeTransformation",
                "CellParameterisation001",
                FatalException,
                ex);
  }
  else
  {
    G4ThreeVector
      origin(fMapCell[copyNo].x(), fMapCell[copyNo].y(), fMapCell[copyNo].z());

    physVol->SetTranslation(origin);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material *
CellParameterisation::ComputeMaterial(const G4int copyNo,
                                      G4VPhysicalVolume *physVol,
                                      const G4VTouchable *)
{
  if (fMaterial[copyNo] == 3) // fMaterial 3 is blue
  {
    physVol->SetName("physicalMat3");
    physVol->GetLogicalVolume()->SetVisAttributes(fBlueAttributes);
    return fBlueMaterial;
  }
  else if (fMaterial[copyNo] == 2) // fMaterial 2 is green
  {
    physVol->SetName("physicalMat2");
    physVol->GetLogicalVolume()->SetVisAttributes(fGreenAttributes);
    return fGreenMaterial;
  }
  else if (fMaterial[copyNo] == 1) // fMaterial 1 is red
  {
    physVol->SetName("physicalMat1");
    physVol->GetLogicalVolume()->SetVisAttributes(fRedAttributes);
    return fRedMaterial;
  }

  return physVol->GetLogicalVolume()->GetMaterial();
}
