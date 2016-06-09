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
#ifndef LXeMainVolume_H
#define LXeMainVolume_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeMainVolume : public G4PVPlacement
{
public:
  LXeMainVolume(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		G4bool pMany,
		G4int pCopyNo,
		LXeDetectorConstruction* c);
private:
  void VisAttributes();
  void SurfaceProperties();

  void PlacePMTs(G4LogicalVolume* pmt_Log,
		 G4RotationMatrix* rot, G4double &a, G4double &b, G4double da,
		 G4double db, G4double amin, G4double bmin, G4int na, G4int nb,
		 G4double &x, G4double &y, G4double &z, G4int &k,LXePMTSD* sd);

  void CopyValues();

  G4bool updated;
  
  LXeDetectorConstruction* constructor;

  G4double scint_x;
  G4double scint_y;
  G4double scint_z;
  G4double d_mtl;
  G4int nx;
  G4int ny;
  G4int nz;
  G4double outerRadius_pmt;
  G4bool sphereOn;
  G4double refl;

  //Basic Volumes
  //
  G4Box* scint_box;
  G4Box* housing_box;
  G4Tubs* pmt;
  G4Tubs* photocath;
  G4Sphere* sphere;


  // Logical volumes
  //
  G4LogicalVolume* scint_log;
  static G4LogicalVolume* housing_log;
  G4LogicalVolume* pmt_log;
  G4LogicalVolume* photocath_log;
  G4LogicalVolume* sphere_log;
  
  // Physical volumes
  //
  G4VPhysicalVolume* scint_phys;
  //keeping pointers to these is pointless really since there are many of them
  //but I'm doing it to be consistent
  G4VPhysicalVolume* pmt_phys;
  G4VPhysicalVolume* photocath_phys;
  G4VPhysicalVolume* sphere_phys;

  //Sensitive Detectors
  static LXeScintSD* scint_SD;
  static LXePMTSD* pmt_SD;
};

#endif
