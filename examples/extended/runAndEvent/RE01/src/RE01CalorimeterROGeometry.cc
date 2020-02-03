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
/// \file runAndEvent/RE01/src/RE01CalorimeterROGeometry.cc
/// \brief Implementation of the RE01CalorimeterROGeometry class
//
<<<<<<< HEAD
// $Id: RE01CalorimeterROGeometry.cc 75295 2013-10-30 09:32:52Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//

#include "RE01CalorimeterROGeometry.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterROGeometry::RE01CalorimeterROGeometry()
  : G4VReadOutGeometry()
{
#include "RE01DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterROGeometry::RE01CalorimeterROGeometry(G4String aString)
  : G4VReadOutGeometry(aString)
{
#include "RE01DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterROGeometry::~RE01CalorimeterROGeometry()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* RE01CalorimeterROGeometry::Build()
{
  //  Material Information imported from NIST database.
  //
  G4NistManager* NISTman = G4NistManager::Instance();

  // A dummy material is used to fill the volumes of the readout geometry.
  // ( It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)
  G4Material* dummyMat  = NISTman->FindOrBuildMaterial("G4_AIR");

  //Builds the ReadOut World:
  G4Box *ROWorldBox = new G4Box("ROWorldBox",fExpHall_x,fExpHall_y,fExpHall_z);
  G4LogicalVolume *ROWorldLog = new G4LogicalVolume(ROWorldBox, dummyMat,
                                                    "ROWorldLogical", 0, 0, 0);
  G4PVPlacement *ROWorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                                 "ROWorldPhysical",
                                                 ROWorldLog,
                                                 0,false,0);
  // Calorimeter volume:
  G4VSolid * caloROtub
    = new G4Tubs("caloROtub",fCaloTubs_rmin,fCaloTubs_rmax,
                 fCaloTubs_dz,fCaloTubs_sphi,fCaloTubs_dphi);
  G4LogicalVolume * caloROlog
    = new G4LogicalVolume(caloROtub,dummyMat,"caloROlogical",0,0,0);
  G4VPhysicalVolume * caloROphys
    = new G4PVPlacement(0,G4ThreeVector(),"calROphysical",caloROlog,
                        ROWorldPhys,false,0);

  // -------------------------------
  // Calorimeter readout division:
  // -------------------------------
  // Phi division first: 48 sectors
  G4VSolid * caloROphiDivisionTub
    = new G4Tubs("caloROphiDivision", fCaloCell_rmin, fCaloCell_rmax,
                 fCaloCell_dz, fCaloCell_sphi, fCaloCell_dphi);
  G4LogicalVolume * caloROphiDivisionLog
    = new G4LogicalVolume(caloROphiDivisionTub, dummyMat, 
                          "caloROphiDivisionLogical",0,0,0);
  G4VPhysicalVolume * caloROphiDivisionPhys
    = new G4PVReplica("caloROphiDivisionPhysical", caloROphiDivisionLog, 
                      caloROphys, kPhi, fSegmentsinPhi, fCaloCell_dphi);
  // then z division: 20 slices:
  G4VSolid * caloROcellTub
    = new G4Tubs("caloROcellTub", fCaloRing_rmin, fCaloRing_rmax,
                 fCaloRing_dz, fCaloRing_sphi, fCaloRing_dphi);
  G4LogicalVolume * caloROcellLog
    = new G4LogicalVolume(caloROcellTub, dummyMat, "caloROcellLogical",0,0,0);
//  G4VPhysicalVolume * caloROcellPhys =
  new G4PVReplica("caloROcellPhysical", caloROcellLog, caloROphiDivisionPhys,
                      kZAxis, fSegmentsinZ,2.*fCaloRing_dz);

  return ROWorldPhys;
}
