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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "CellParameterisation.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Region.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:

    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume *Construct() override;

    inline auto *GetLogicalMedium() const { return fLogicMedium; };

    void SetTargetMaterial(const G4String&);

    void SetRedDensity(const G4double&);
    void SetGreenDensity(const G4double&);
    void SetBlueDensity(const G4double&);

    void SetShiftX(const G4double&);
    void SetShiftY(const G4double&);
    void SetShiftZ(const G4double&);

    void SetMediumSizeXY(const G4double&);
    void SetMediumSizeZ(const G4double&);

    void SetWorldSizeXY(const G4double&);
    void SetWorldSizeZ(const G4double&);

    void SetPhantomFileName(const G4String&);

  private:

    void DefineMaterials();

    G4VPhysicalVolume *ConstructLine();

    G4double fDensityRed = 1.0;
    G4double fDensityGreen = 1.0;
    G4double fDensityBlue = 1.0;

    G4double fShiftX = 0.*um;
    G4double fShiftY = 0.*um;
    G4double fShiftZ = 0.*um;

    G4double fWorldSizeXY = 0.;
    G4double fWorldSizeZ = 0.;

    G4double fMediumSizeXY = 0.;
    G4double fMediumSizeZ = 0.;

    G4Material *fDefaultMaterial = nullptr;
    G4Material *fMediumMaterial = nullptr;
    G4Material *fRedMaterial = nullptr;
    G4Material *fGreenMaterial = nullptr;
    G4Material *fBlueMaterial = nullptr;
    G4Material *fPhantomMaterial = nullptr;

    G4VPhysicalVolume *fPhysiWorld = nullptr;
    G4LogicalVolume *fLogicWorld = nullptr;
    G4Box *fSolidWorld = nullptr;

    G4VPhysicalVolume *fPhysiMedium = nullptr;
    G4LogicalVolume *fLogicMedium = nullptr;
    G4Box *fSolidMedium = nullptr;

    G4VPhysicalVolume *fPhysiPhantom = nullptr;
    G4LogicalVolume *fLogicPhantom = nullptr;
    G4Box *fSolidPhantom = nullptr;
    CellParameterisation *fPhantomParam = nullptr;

    DetectorMessenger* fDetectorMessenger = nullptr;

    G4String fPhantomFileName = "";
    G4Region* fPhantomRegion = nullptr;
};

#endif
