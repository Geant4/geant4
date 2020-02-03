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
/// \file eventgenerator/HepMC/HepMCEx01/src/ExN04CalorimeterROGeometry.cc
/// \brief Implementation of the ExN04CalorimeterROGeometry class
//
//

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "ExN04CalorimeterROGeometry.hh"
#include "ExN04DummySD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04CalorimeterROGeometry::ExN04CalorimeterROGeometry()
  : G4VReadOutGeometry()
{
#include "ExN04DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04CalorimeterROGeometry::ExN04CalorimeterROGeometry(G4String aString)
  : G4VReadOutGeometry(aString)
{
#include "ExN04DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN04CalorimeterROGeometry::~ExN04CalorimeterROGeometry()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExN04CalorimeterROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // ( It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)
  G4Material* dummyMat  =
     new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  //Builds the ReadOut World:
  G4Box* ROWorldBox = new G4Box("ROWorldBox",
                                fexpHall_x, fexpHall_y, fexpHall_z);
  G4LogicalVolume* ROWorldLog = new G4LogicalVolume(ROWorldBox, dummyMat,
                                                    "ROWorldLogical", 0, 0, 0);
  G4PVPlacement* ROWorldPhys = new G4PVPlacement(0, G4ThreeVector(),
                                                 "ROWorldPhysical",
                                                 ROWorldLog,
                                                 0, false, 0);
  // Calorimeter volume:
  G4VSolid* caloROtub
    = new G4Tubs("caloROtub", fcaloTubs_rmin, fcaloTubs_rmax,
                 fcaloTubs_dz, fcaloTubs_sphi, fcaloTubs_dphi);
  G4LogicalVolume* caloROlog
    = new G4LogicalVolume(caloROtub, dummyMat, "caloROlogical",0,0,0);
  G4VPhysicalVolume* caloROphys
    = new G4PVPlacement(0, G4ThreeVector(), "calROphysical", caloROlog,
                        ROWorldPhys, false, 0);

  // -------------------------------
  // Calorimeter readout division:
  // -------------------------------
  // Phi division first: 48 sectors
  G4VSolid* caloROphiDivisionTub
    = new G4Tubs("caloROphiDivision", fcaloCell_rmin, fcaloCell_rmax,
                 fcaloCell_dz, fcaloCell_sphi, fcaloCell_dphi);
  G4LogicalVolume* caloROphiDivisionLog
    = new G4LogicalVolume(caloROphiDivisionTub,
                          dummyMat, "caloROphiDivisionLogical",0,0,0);
  G4VPhysicalVolume* caloROphiDivisionPhys
    = new G4PVReplica("caloROphiDivisionPhysical", caloROphiDivisionLog,
                      caloROphys, kPhi, fsegmentsinPhi, fcaloCell_dphi);
  // then z division: 20 slices:
  G4VSolid* caloROcellTub
    = new G4Tubs("caloROcellTub", fcaloRing_rmin, fcaloRing_rmax,
                 fcaloRing_dz, fcaloRing_sphi, fcaloRing_dphi);
  G4LogicalVolume* caloROcellLog
    = new G4LogicalVolume(caloROcellTub, dummyMat, "caloROcellLogical",0,0,0);
  //  G4VPhysicalVolume * caloROcellPhys =
  new G4PVReplica("caloROcellPhysical", caloROcellLog, caloROphiDivisionPhys,
                  kZAxis, fsegmentsinZ, 2.*fcaloRing_dz);


  //Flags the cells as sensitive .The pointer here serves
  // as a flag only to check for sensitivity.
  // (Could we make it by a simple cast of a non-NULL value ?)
  ExN04DummySD * dummySensi = new ExN04DummySD;
  caloROcellLog->SetSensitiveDetector(dummySensi);

  return ROWorldPhys;
}
