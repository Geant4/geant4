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
/// \file runAndEvent/RE04/src/RE04ParallelWorldConstruction.cc
/// \brief Implementation of the RE04ParallelWorldConstruction class
//
//
#include "RE04ParallelWorldConstruction.hh"
#include "RE04ParallelWorldParam.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ParallelWorldConstruction
::RE04ParallelWorldConstruction(G4String& parallelWorldName)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ParallelWorldConstruction::~RE04ParallelWorldConstruction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04ParallelWorldConstruction::Construct()
{
  if(fConstructed) return;
  fConstructed = true;

  //
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

  //
  // material defined in the mass world
  //
  G4Material* water = G4Material::GetMaterial("G4_WATER");

  //
  // parallel world placement box
  //
  G4VSolid* paraBox = new G4Box("paraBox",5.0*cm,30.0*cm,5.0*cm);
  G4LogicalVolume* paraBoxLogical = new G4LogicalVolume(paraBox,water,
                                                        "paraBox");
  new G4PVPlacement(0,G4ThreeVector(-25.0*cm,0.,0.),paraBoxLogical,
                    "paraBox",worldLogical,false,0);

  //
  // mother of parallel world parameterized volumes
  //
  G4VSolid* paraMom = new G4Box("paraMom",20.0*cm,40.0*cm,20.0*cm);
  G4LogicalVolume* paraMomLogical = new G4LogicalVolume(paraMom,0,"paraMom");
  new G4PVPlacement(0,G4ThreeVector(10.0*cm,0.,0.),paraMomLogical,"paraMom",
                    worldLogical,false,0);

  //
  // parallel world parameterized volumes
  //
  G4VSolid* paraPara = new G4Box("paraPara",5.0*cm,15.0*cm,10.0*cm);
  G4LogicalVolume* paraParaLogical = new G4LogicalVolume(paraPara,water,
                                                         "paraPara");
  RE04ParallelWorldParam* param = new RE04ParallelWorldParam();
  new G4PVParameterised("paraPara",paraParaLogical,paraMomLogical,
                           kXAxis, 2, param);

}
