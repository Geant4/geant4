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
/// \file exoticphysics/phonon/src/XDetectorConstruction.cc
/// \brief Implementation of the XDetectorConstruction class
//
// $Id$
//

#include "XDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include "XAlminumElectrodeSensitivity.hh"
#include "XPhysicalLattice.hh"
#include "XLogicalLattice.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

const G4String XDetectorConstruction::fsCrystalMapsDir = 
     getenv("CRYSTALMAPS") ? getenv("CRYSTALMAPS") : "./CrystalMaps";


XDetectorConstruction::XDetectorConstruction():fConstructed(false),fIfField(true)
{
  fLiquidHelium = NULL;
  fGermanium = NULL;
  fAlminum = NULL;
  fTungsten = NULL;
  fWorldPhys = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XDetectorConstruction::~XDetectorConstruction()
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XDetectorConstruction::Construct()
{
  if(!fConstructed)
  { 
    fConstructed = true;
    DefineMaterials();
    SetupGeometry();
  }
  return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XDetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected.......
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fAlminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XDetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
  worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",0,false,0);
  
  //                               
  // Germanium cylinder - this is the volume in which we will propagate phonons
  //  
  G4VSolid* fGermaniumSolid = new G4Tubs("fGermaniumSolid",0.*cm,3.81*cm,1.27*cm, 0.*deg, 360.*deg);
  G4LogicalVolume* fGermaniumLogical = new G4LogicalVolume(fGermaniumSolid,fGermanium,"fGermaniumLogical");
  G4VPhysicalVolume* GePhys = new G4PVPlacement(0,G4ThreeVector(),fGermaniumLogical,"fGermaniumPhysical",worldLogical,false,0);

  //
  //Germanium lattice information
  //
  
  //Get pointer to new logical lattice object
  XLogicalLattice* GeLogical = new XLogicalLattice();
  
  //Load maps for mapping phonon momentum direciton onto propagation direction
  //Convention for polarization state: 0=LON, 1=ST, 2=FT
  if(GeLogical->Load_NMap(161, 321, 0, fsCrystalMapsDir + "/LVec.ssv")){
    if(GeLogical->Load_NMap(161, 321, 1,fsCrystalMapsDir + "/STVec.ssv")){
      if(GeLogical->Load_NMap(161, 321, 2,fsCrystalMapsDir + "/FTVec.ssv")){
        G4cout<<"\nXDetectorConstruction::Loaded all three maps";}}}

  //Load maps for mapping phonon momentum direction onto velocity
  //Convention for polarization state: 0=LON, 1=ST, 2=FT
  if(GeLogical->LoadMap(161, 321, 0, fsCrystalMapsDir +"/L.ssv")){
    if(GeLogical->LoadMap(161, 321, 1, fsCrystalMapsDir +"/ST.ssv")){
      if(GeLogical->LoadMap(161, 321, 2, fsCrystalMapsDir +"/FT.ssv")){
        G4cout<<"\nXDetectorConstruction::Loaded all three velocity maps";}}}
  
  //Set Ge lattice dynamical information
  GeLogical->SetDynamicalConstants(-0.732, -0.708, 0.376, 0.561);
  GeLogical->SetScatteringConstant(3.67e-41*s*s*s);
  GeLogical->SetAnhDecConstant(1.6456e-54*s*s*s*s);
  GeLogical->SetLDOS(0.097834);
  GeLogical->SetSTDOS(0.53539);
  GeLogical->SetFTDOS(0.36677);

  //XPhysicalLattice associates XLogicalLattice with a PhysicalVolume.
  //XLatticeManager3 gives physics processes access to lattices
  XPhysicalLattice* GePhysical = new XPhysicalLattice(GePhys,GeLogical);
  XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
  LM->RegisterLattice(GePhysical);

  //
  // Alminum - crystal end caps. This is where phonon hits are registered
  //
  G4VSolid* fAlminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm, 0.*deg, 360.*deg);

  G4LogicalVolume* fAlminumLogical = new G4LogicalVolume(fAlminumSolid,fAlminum,"fAlminumLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,1.28*cm),fAlminumLogical,"fAlminumPhysical",worldLogical,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-1.28*cm),fAlminumLogical,"fAlminumPhysical",worldLogical,false,1);


  //
  // detector -- Note : Alminum electrode sensitivity is attached to Germanium 
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  XAlminumElectrodeSensitivity* electrodeSensitivity = new XAlminumElectrodeSensitivity("XAlminumElectrode");
  SDman->AddNewDetector(electrodeSensitivity);
  fGermaniumLogical->SetSensitiveDetector(electrodeSensitivity);

  //                                        
  // Visualization attributes
  //
  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  fGermaniumLogical->SetVisAttributes(simpleBoxVisAtt);
  fAlminumLogical->SetVisAttributes(simpleBoxVisAtt);

}


