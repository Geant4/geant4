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
/// \file optical/LXe/include/LXeDetectorConstruction.hh
/// \brief Definition of the LXeDetectorConstruction class
//
//
#ifndef LXeDetectorConstruction_h
#define LXeDetectorConstruction_h 1

#include "LXeDetectorMessenger.hh"

#include "G4Cache.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>

class LXeMainVolume;
class LXePMTSD;
class LXeScintSD;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  LXeDetectorConstruction();
  ~LXeDetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  // Functions to modify the geometry
  void SetDimensions(G4ThreeVector);
  void SetHousingThickness(G4double);
  void SetNX(G4int);
  void SetNY(G4int);
  void SetNZ(G4int);
  void SetPMTRadius(G4double);
  void SetDefaults();
  void SetSaveThreshold(G4int);

  // Get values
  G4int GetNX() const { return fNx; };
  G4int GetNY() const { return fNy; };
  G4int GetNZ() const { return fNz; };
  G4int GetSaveThreshold() const { return fSaveThreshold; };
  G4double GetScintX() const { return fScint_x; }
  G4double GetScintY() const { return fScint_y; }
  G4double GetScintZ() const { return fScint_z; }
  G4double GetHousingThickness() const { return fD_mtl; }
  G4double GetPMTRadius() const { return fOuterRadius_pmt; }
  G4double GetSlabZ() const { return fSlab_z; }

  void SetSphereOn(G4bool);
  static G4bool GetSphereOn() { return fSphereOn; }

  void SetHousingReflectivity(G4double);
  G4double GetHousingReflectivity() const { return fRefl; }

  void SetWLSSlabOn(G4bool b);
  G4bool GetWLSSlabOn() const { return fWLSslab; }

  void SetMainVolumeOn(G4bool b);
  G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

  void SetNFibers(G4int n);
  G4int GetNFibers() const { return fNfibers; }

  void SetMainScintYield(G4double);
  void SetWLSScintYield(G4double);

 private:
  void DefineMaterials();

  LXeDetectorMessenger* fDetectorMessenger = nullptr;

  G4Box* fExperimentalHall_box = nullptr;
  G4LogicalVolume* fExperimentalHall_log = nullptr;
  G4VPhysicalVolume* fExperimentalHall_phys = nullptr;

  // Materials & Elements
  G4Element* fN = nullptr;
  G4Element* fO = nullptr;
  G4Element* fC = nullptr;
  G4Element* fH = nullptr;
  G4Material* fLXe = nullptr;
  G4Material* fAl = nullptr;
  G4Material* fAir = nullptr;
  G4Material* fVacuum = nullptr;
  G4Material* fGlass = nullptr;
  G4Material* fPstyrene = nullptr;
  G4Material* fPMMA = nullptr;
  G4Material* fPethylene1 = nullptr;
  G4Material* fPethylene2 = nullptr;

  // Geometry
  G4double fScint_x = 17.8 * CLHEP::cm;
  G4double fScint_y = 17.8 * CLHEP::cm;
  G4double fScint_z = 22.6 * CLHEP::cm;
  G4double fD_mtl = 0.0635 * CLHEP::cm;
  G4int fNx = 2;
  G4int fNy = 2;
  G4int fNz = 3;
  G4int fSaveThreshold = 0;
  G4double fOuterRadius_pmt = 2.3 * CLHEP::cm;
  G4int fNfibers = 15;
  static G4bool fSphereOn;
  G4double fRefl = 1.;
  G4bool fWLSslab = false;
  G4bool fMainVolumeOn = true;
  G4double fSlab_z = 2.5 * CLHEP::mm;

  LXeMainVolume* fMainVolume = nullptr;

  G4MaterialPropertiesTable* fLXe_mt = nullptr;
  G4MaterialPropertiesTable* fMPTPStyrene = nullptr;

  // Sensitive Detectors
  G4Cache<LXeScintSD*> fScint_SD;
  G4Cache<LXePMTSD*> fPmt_SD;
};

#endif
