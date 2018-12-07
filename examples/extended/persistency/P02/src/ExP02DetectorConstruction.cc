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
/// \file persistency/P02/src/ExP02DetectorConstruction.cc
/// \brief Implementation of the ExP02DetectorConstruction class
//
//
//
#include "ExP02DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"

#include "ExP02GeoTree.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP02DetectorConstruction::ExP02DetectorConstruction()
  :  G4VUserDetectorConstruction(),
    fExperimentalHall_log(0), fTracker_log(0),
    fCalorimeterBlock_log(0), fCalorimeterLayer_log(0),
    fExperimentalHall_phys(0), fCalorimeterLayer_phys(0),
    fCalorimeterBlock_phys(0), fTracker_phys(0)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP02DetectorConstruction::~ExP02DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExP02DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4Material* Ar = 
  new G4Material("ArgonGas", z= 18., a= 39.95*g/mole, density= 1.782*mg/cm3);

  G4Material* Al = 
  new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  G4Material* Pb = 
  new G4Material("Lead", z= 82., a= 207.19*g/mole, density= 11.35*g/cm3);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  //------------------------------ beam line along x axis

  G4double expHall_x = 3.0*m;
  G4double expHall_y = 1.0*m;
  G4double expHall_z = 1.0*m;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Ar,"expHall_log",0,0,0);
  fExperimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                                      fExperimentalHall_log,"expHall",0,false,0);

  //------------------------------ a tracker tube

  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 60.*cm;
  G4double hightOfTheTube = 50.*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
                                    outerRadiusOfTheTube,hightOfTheTube,
                                    startAngleOfTheTube,spanningAngleOfTheTube);
  fTracker_log = new G4LogicalVolume(tracker_tube,Al,"tracker_log",0,0,0);
  G4double trackerPos_x = -1.0*m;
  G4double trackerPos_y = 0.*m;
  G4double trackerPos_z = 0.*m;
  fTracker_phys = new G4PVPlacement(0,
             G4ThreeVector(trackerPos_x,trackerPos_y,trackerPos_z),
             fTracker_log,"tracker",fExperimentalHall_log,false,0);

  //------------------------------ a calorimeter block

  G4double block_x = 1.0*m;
  G4double block_y = 50.0*cm;
  G4double block_z = 50.0*cm;
  G4Box* calorimeterBlock_box = new G4Box("calBlock_box",block_x,
                                          block_y,block_z);
  fCalorimeterBlock_log = new G4LogicalVolume(calorimeterBlock_box,
                                             Pb,"caloBlock_log",0,0,0);
  G4double blockPos_x = 1.0*m;
  G4double blockPos_y = 0.0*m;
  G4double blockPos_z = 0.0*m;
  fCalorimeterBlock_phys = new G4PVPlacement(0,
             G4ThreeVector(blockPos_x,blockPos_y,blockPos_z),
             fCalorimeterBlock_log,"caloBlock",fExperimentalHall_log,false,0);

  //------------------------------ calorimeter layers

  G4double calo_x = 1.*cm;
  G4double calo_y = 40.*cm;
  G4double calo_z = 40.*cm;
  G4Box* calorimeterLayer_box = new G4Box("caloLayer_box",
                                          calo_x,calo_y,calo_z);
  fCalorimeterLayer_log = new G4LogicalVolume(calorimeterLayer_box,
                                             Al,"caloLayer_log",0,0,0);
  for(G4int i=0;i<19;i++) // loop for 19 layers
  {
    G4double caloPos_x = (i-9)*10.*cm;
    G4double caloPos_y = 0.0*m;
    G4double caloPos_z = 0.0*m;
    fCalorimeterLayer_phys = new G4PVPlacement(0,
               G4ThreeVector(caloPos_x,caloPos_y,caloPos_z),
               fCalorimeterLayer_log,"caloLayer",fCalorimeterBlock_log,false,i);
  }

  //------------------------------------------------------------------

  // writing the geometry to a ROOT file

  // initialize ROOT
  TSystem ts;

  gSystem->Load("libExP02ClassesDict");
  
  //  gDebug = 1;

  const G4ElementTable* eltab = G4Element::GetElementTable();

  ExP02GeoTree* geotree = new ExP02GeoTree(fExperimentalHall_phys, eltab);

  TFile fo("geo.root","RECREATE");

  fo.WriteObject(geotree, "my_geo");
  
  return fExperimentalHall_phys;
}

