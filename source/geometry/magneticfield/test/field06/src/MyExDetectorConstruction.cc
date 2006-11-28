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
// $Id: MyExDetectorConstruction.cc,v 1.1 2006-11-28 15:19:42 tnikitin Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

ExN01DetectorConstruction::ExN01DetectorConstruction()
 :  experimentalHall_log(0), tracker_log(0),
    calorimeterBlock_log(0), calorimeterLayer_log(0),
    experimentalHall_phys(0), calorimeterLayer_phys(0),
    calorimeterBlock_phys(0), tracker_phys(0)
{;}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4double pressure;
  G4double temperature;
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4Material* Ar = 
  new G4Material("ArgonGas", z= 18., a= 39.95*g/mole, density= 1.782*mg/cm3);

  G4Material* Al = 
  new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  G4Material* Pb = 
  new G4Material("Lead", z= 82., a= 207.19*g/mole, density= 11.35*g/cm3);
//
// What about vacuum ?  Vacuum is an ordinary gas with very low density
//

density     = universe_mean_density;                //from PhysicalConstants.h
pressure    = 1.e-19*pascal;
temperature = 0.1*kelvin;
G4Material* Galactic=new G4Material("Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

//
// define a material from elements.   case 2: mixture by fractional mass
//
 //This function illustrates the possible ways to define materials
 
G4String symbol;         




////define Elements



G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);

G4Material* Air = 
new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);



  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)
  //------------------------------ beam line along x axis

  G4double expHall_x = 200.0*cm;
  G4double expHall_y = 200.0*cm;
  G4double expHall_z = 200.0*cm;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             Galactic,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                                      experimentalHall_log,"expHall",0,false,0);

  //------------------------------ a tracker tube

  G4double innerRadiusOfTheTube = 0.*cm;
  G4double outerRadiusOfTheTube = 70.*cm;
  G4double hightOfTheTube = 70.*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;
  G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
                                    outerRadiusOfTheTube,hightOfTheTube,
                                    startAngleOfTheTube,spanningAngleOfTheTube);
  tracker_log = new G4LogicalVolume(tracker_tube,Ar,"tracker_log",0,0,0);
  G4double trackerPos_x = -0.0*cm;
  G4double trackerPos_y = 0.*cm;
  G4double trackerPos_z = 0.*cm;
  tracker_phys = new G4PVPlacement(0,
             G4ThreeVector(trackerPos_x,trackerPos_y,trackerPos_z),
             tracker_log,"tracker",experimentalHall_log,false,0);


  //------------------------------ a calorimeter block

  // G4double block_x = 1.0*m;
  //G4double block_y = 50.0*cm;
  //G4double block_z = 50.0*cm;
  //G4Box* calorimeterBlock_box = new G4Box("calBlock_box",block_x,
  //                                        block_y,block_z);
  //calorimeterBlock_log = new G4LogicalVolume(calorimeterBlock_box,
  //                                           Pb,"caloBlock_log",0,0,0);
  //G4double blockPos_x = 1.0*m;
  //G4double blockPos_y = 0.0*m;
  //G4double blockPos_z = 0.0*m;
  //calorimeterBlock_phys = new G4PVPlacement(0,
  //           G4ThreeVector(blockPos_x,blockPos_y,blockPos_z),
  //           calorimeterBlock_log,"caloBlock",experimentalHall_log,false,0);

  //------------------------------ calorimeter layers

  //G4double calo_x = 1.*cm;
  //G4double calo_y = 40.*cm;
  //G4double calo_z = 40.*cm;
  //G4Box* calorimeterLayer_box = new G4Box("caloLayer_box",
  //                                        calo_x,calo_y,calo_z);
  //calorimeterLayer_log = new G4LogicalVolume(calorimeterLayer_box,
  //                                           Al,"caloLayer_log",0,0,0);
  //for(G4int i=0;i<19;i++) // loop for 19 layers
  //{
  //  G4double caloPos_x = (i-9)*10.*cm;
  //  G4double caloPos_y = 0.0*m;
  //  G4double caloPos_z = 0.0*m;
  //  calorimeterLayer_phys = new G4PVPlacement(0,
  //             G4ThreeVector(caloPos_x,caloPos_y,caloPos_z),
  //             calorimeterLayer_log,"caloLayer",calorimeterBlock_log,false,i);
  //}

  //------------------------------------------------------------------

  return experimentalHall_phys;
}

