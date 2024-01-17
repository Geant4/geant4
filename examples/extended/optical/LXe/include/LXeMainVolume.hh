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
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef LXeMainVolume_h
#define LXeMainVolume_h 1

#include "LXeDetectorConstruction.hh"

#include "G4PVPlacement.hh"

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class LXeMainVolume : public G4PVPlacement
{
 public:
  LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                LXeDetectorConstruction* c);

  G4LogicalVolume* GetLogPhotoCath() { return fPhotocath_log; }
  G4LogicalVolume* GetLogScint() { return fScint_log; }

  std::vector<G4ThreeVector> GetPmtPositions() { return fPmtPositions; }

 private:
  void VisAttributes();
  void SurfaceProperties();

  void PlacePMTs(G4LogicalVolume* pmt_Log, G4RotationMatrix* rot, G4double& a,
                 G4double& b, G4double da, G4double db, G4double amin,
                 G4double bmin, G4int na, G4int nb, G4double& x, G4double& y,
                 G4double& z, G4int& k);

  void CopyValues();

  LXeDetectorConstruction* fConstructor = nullptr;

  G4double fScint_x = 0.;
  G4double fScint_y = 0.;
  G4double fScint_z = 0.;
  G4double fD_mtl = 0.;
  G4int fNx = 0;
  G4int fNy = 0;
  G4int fNz = 0;
  G4double fOuterRadius_pmt = 0.;
  G4bool fSphereOn = false;
  G4double fRefl = 0.;

  // Basic Volumes
  //
  G4Box* fScint_box = nullptr;
  G4Box* fHousing_box = nullptr;
  G4Tubs* fPmt = nullptr;
  G4Tubs* fPhotocath = nullptr;
  G4Sphere* fSphere = nullptr;

  // Logical volumes
  //
  G4LogicalVolume* fScint_log = nullptr;
  G4LogicalVolume* fHousing_log = nullptr;
  G4LogicalVolume* fPmt_log = nullptr;
  G4LogicalVolume* fPhotocath_log = nullptr;
  G4LogicalVolume* fSphere_log = nullptr;

  // Sensitive Detectors positions
  std::vector<G4ThreeVector> fPmtPositions;
};

#endif
