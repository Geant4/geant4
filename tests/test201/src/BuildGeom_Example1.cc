// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildGeom_Example1.cc,v 1.2 1999-12-15 14:55:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Makoto Asai - based on Long Baseline Neutrino Observatory experiment.
// "Complete" detector but no rotated volumes.

#include "BuildGeom_Example1.hh"

#include <math.h>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Navigator.hh"
#include "G4GeometryManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"


G4VPhysicalVolume* BuildGeom_Example1()
{

//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", iz=14., a);
  a = 207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  density = 5.2*g/cm3;
  G4Material* LeadGlass = new G4Material(name="LeadGlass", density, nel=3);
  LeadGlass->AddElement(elO, .199);
  LeadGlass->AddElement(elSi, .127);
  LeadGlass->AddElement(elPb, .674);
  density = 1.0*g/cm3;
  G4Material* Water = new G4Material(name="water",density,nel=2);
  Water->AddElement(elH, (G4int)2);
  Water->AddElement(elO, (G4int)1);
  density = 1.032*g/cm3;
  G4Material* Scinti = new G4Material(name="Scinti",density,nel=2);
  Scinti->AddElement(elH, (G4int)1);
  Scinti->AddElement(elC, (G4int)1);


// Experimental hall (world volume)

  G4Box *myWorldBox= new G4Box("WBox",600.*cm,1200.*cm,1000.*cm);
  G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,Air,
						    "WLog", 0, 0, 0);
  myWorldLog -> SetVisAttributes (&G4VisAttributes::Invisible);
  G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);


  G4double waterSciFiYpos = 100.*cm;

  G4VisAttributes* red   = new G4VisAttributes (G4Colour (1,0,0));
  G4VisAttributes* green = new G4VisAttributes (G4Colour (0,1,0));
  G4VisAttributes* blue  = new G4VisAttributes (G4Colour (0,0,1));

// WaterSciFi overall

  G4Box* watSFBox = new G4Box("WSBox",160.*cm,65.*cm,160.*cm);
  G4LogicalVolume* watSFLog = new G4LogicalVolume(watSFBox,Al,"WSLog",0,0,0);
  watSFLog -> SetVisAttributes (red);
  G4PVPlacement* watSFPhys = new G4PVPlacement(0,
    G4ThreeVector(0.*cm,waterSciFiYpos,0.*cm),
    "watSFPhys",watSFLog,myWorldPhys,false,0);

// WaterSciFi water box

  G4Box* watSFWBox = new G4Box("WSWBox",159.*cm,64.*cm,159.*cm);
  G4LogicalVolume* watSFWLog = 
    new G4LogicalVolume(watSFWBox,Water,"WSWLog",0,0,0);
  watSFWLog -> SetVisAttributes (green);
  G4PVPlacement* watSFWPhys = new G4PVPlacement(0,
	  G4ThreeVector(0.*cm,0.*cm,0.*cm),
	  "watSFWPhys",watSFWLog,watSFPhys,false,0);

// Scinti plate (XXYY are merged into one plane)

  G4Box* watSFSBox = new G4Box("WSSBox",159.*cm,0.1*cm,159.*cm);
  G4LogicalVolume* watSFSLog = new G4LogicalVolume(watSFSBox,Scinti,"WSSLog",0,0,0);
  watSFSLog -> SetVisAttributes (blue);
  G4PVPlacement* watSFSPhys;
  for(int iSFS=-10; iSFS<=10; iSFS++)
  {
    G4double SFSypos = iSFS*6.0*cm;
    watSFSPhys = new G4PVPlacement(0,
	G4ThreeVector(0.*cm,SFSypos,0.*cm),
	"watSFSPhys",watSFSLog,watSFWPhys,false,iSFS+10);
  }

// LeadGlass support unit

  G4Tubs* calUnitTubsR 
    = new G4Tubs("CUTubsR",174.636*cm,213.9*cm,30.5*cm,45.0*deg,45.0*deg);
  G4LogicalVolume* calUnitLogR
    = new G4LogicalVolume(calUnitTubsR,Fe,"CUlogR",0,0,0);
  G4Tubs* calUnitTubsL 
    = new G4Tubs("CUTubsL",174.636*cm,213.9*cm,30.5*cm,90.0*deg,45.0*deg);
  G4LogicalVolume* calUnitLogL
    = new G4LogicalVolume(calUnitTubsL,Fe,"CUlogL",0,0,0);

  G4PVPlacement* calUnitPhysR;
  G4PVPlacement* calUnitPhysL;
  for(int iCalUnit=-2; iCalUnit<=2; iCalUnit++)
  {
    G4double zCalRow = iCalUnit*61.0*cm;
    calUnitPhysR = new G4PVPlacement(0,
               G4ThreeVector(0.*cm,waterSciFiYpos,zCalRow),
               "calUnitPhysR",calUnitLogR,myWorldPhys,false,iCalUnit+2);
    calUnitPhysL = new G4PVPlacement(0,
               G4ThreeVector(0.*cm,waterSciFiYpos,zCalRow),
               "calUnitPhysL",calUnitLogL,myWorldPhys,false,iCalUnit+2);
  }

// LeadGlass body

  G4Tubs* calRowTubsR 
    = new G4Tubs("CRTubsR",174.636*cm,212.9*cm,30.5*cm,45.0*deg,45.0*deg);
  G4LogicalVolume* calRowLogR
    = new G4LogicalVolume(calRowTubsR,LeadGlass,"CRlogR",0,0,0);
  G4PVPlacement* calRowPhysR = new G4PVPlacement(0,
	       G4ThreeVector(0.*cm,0.*cm,0.*cm),
	       "calRowPhysR",calRowLogR,calUnitPhysR,false,0);
  G4Tubs* calRowTubsL 
    = new G4Tubs("CRTubsL",174.636*cm,212.9*cm,30.5*cm,90.0*deg,45.0*deg);
  G4LogicalVolume* calRowLogL
    = new G4LogicalVolume(calRowTubsL,LeadGlass,"CRlogL",0,0,0);
  G4PVPlacement* calRowPhysL = new G4PVPlacement(0,
	       G4ThreeVector(0.*cm,0.*cm,0.*cm),
	       "calRowPhysL",calRowLogL,calUnitPhysL,false,0);

// Muon filter Iron

  G4Box* muFeBox = new G4Box("muFe",400.*cm,5.0*cm,400.*cm);
  G4LogicalVolume* muFeLog = new G4LogicalVolume(muFeBox,Fe,"MFLog",0,0,0);
  G4Box* muChBox = new G4Box("muCh",400.*cm,10.0*cm,400.*cm);
  G4LogicalVolume* muChLog = new G4LogicalVolume(muChBox,Air,"MCLog",0,0,0);

  G4double muonChamberYPos = 300.*cm+waterSciFiYpos;

  G4int muFeNo = 0;
  G4int muChNo = 0;
  G4PVPlacement* muFePhys;
  G4PVPlacement* muChPhys;

  for(G4int iMC1=0; iMC1<4; iMC1++)
  {
    muonChamberYPos += 10.*cm;
    muChPhys = new G4PVPlacement(0,
     G4ThreeVector(0.*cm,muonChamberYPos,0.*cm),
     "muChPhys",muChLog,myWorldPhys,false,muChNo++);
    muonChamberYPos += 10.*cm;

    muonChamberYPos += 5.*cm;
    muFePhys = new G4PVPlacement(0,
     G4ThreeVector(0.*cm,muonChamberYPos,0.*cm),
     "muFePhys",muFeLog,myWorldPhys,false,muFeNo++);
    muonChamberYPos += 5.*cm;
  }

  for(G4int iMC2=0; iMC2<8; iMC2++)
  {
    muonChamberYPos += 10.*cm;
    muChPhys = new G4PVPlacement(0,
     G4ThreeVector(0.*cm,muonChamberYPos,0.*cm),
     "muChPhys",muChLog,myWorldPhys,false,muChNo++);
    muonChamberYPos += 10.*cm;

    muonChamberYPos += 5.*cm;
    muFePhys = new G4PVPlacement(0,
     G4ThreeVector(0.*cm,muonChamberYPos,0.*cm),
     "muFePhys",muFeLog,myWorldPhys,false,muFeNo++);
    muonChamberYPos += 5.*cm;

    muonChamberYPos += 5.*cm;
    muFePhys = new G4PVPlacement(0,
     G4ThreeVector(0.*cm,muonChamberYPos,0.*cm),
     "muFePhys",muFeLog,myWorldPhys,false,muFeNo++);
    muonChamberYPos += 5.*cm;
  }

  return myWorldPhys;
}

