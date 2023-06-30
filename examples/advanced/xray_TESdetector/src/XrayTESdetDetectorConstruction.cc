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
/// \file XrayTESdetDetectorConstruction.cc
/// \brief Implementation of the XrayTESdetDetectorConstruction class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetDetectorConstruction.hh"
#include "XrayTESdetDetectorMessenger.hh"
#include "XrayTESdetDetParameterisation.hh"
#include "G4Material.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4PVParameterised.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Hype.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// #include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include "G4GDMLParser.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"


XrayTESdetDetectorConstruction::XrayTESdetDetectorConstruction()
: fStepLimitPolyfilta(nullptr),fStepLimitAlfilta(nullptr),fStepLimitmembrane(nullptr),fExperimentalHall_phys(nullptr),fDetectorMessenger(nullptr)
{
  fReadFile = "";
  fDetectorMessenger = new XrayTESdetDetectorMessenger(this);
}


XrayTESdetDetectorConstruction::~XrayTESdetDetectorConstruction()
{
  delete fStepLimitAlfilta;
  delete fStepLimitPolyfilta;
  delete fStepLimitmembrane;
  delete fDetectorMessenger;
}


G4VPhysicalVolume* XrayTESdetDetectorConstruction::Construct()
{
  return ConstructDetector();
}


G4VPhysicalVolume* XrayTESdetDetectorConstruction::ConstructDetector()
{
  G4VPhysicalVolume* fWorld_phys;
  auto *InnerRegion = new G4Region("InnerRegion");

  if (fReadFile == "")
  {
    G4cout << "Build a new geometry" << G4endl;

    // -------------------------------------------------------------------------------------------
    // --------------- Define materials ----------------------------------------------------------
    // -------------------------------------------------------------------------------------------
    G4String symbol;
    G4int ncomponents;
    G4double density;
    G4double fractionmass;

    auto *man = G4NistManager::Instance();
    auto *Si = man->FindOrBuildMaterial("G4_Si");
    auto *Cu = man->FindOrBuildMaterial("G4_Cu");
    auto *Bi = man->FindOrBuildMaterial("G4_Bi");
    auto *Al = man->FindOrBuildMaterial("G4_Al");
    auto *Ni = man->FindOrBuildMaterial("G4_Ni");
    auto *Au = man->FindOrBuildMaterial("G4_Au");
    auto *N = man->FindOrBuildMaterial("G4_N");
    auto *Fe = man->FindOrBuildMaterial("G4_Fe");
    auto *C = man->FindOrBuildMaterial("G4_C");
    auto *Mn = man->FindOrBuildMaterial("G4_Mn");
    auto *Ph = man->FindOrBuildMaterial("G4_P");
    auto *S = man->FindOrBuildMaterial("G4_S");
    auto *Cr = man->FindOrBuildMaterial("G4_Cr");
    auto *Mo = man->FindOrBuildMaterial("G4_Mo");
    auto *Ti = man->FindOrBuildMaterial("G4_Ti");

    // Define the material of the ring
    auto *RingAl = new G4Material("RingAl", density=0.736*g/cm3, ncomponents=1);
    RingAl->AddMaterial(Al, fractionmass=1);

    // Define cryoperm
    auto *Cryoperm = new G4Material("Cryoperm", density=8.7*g/cm3, ncomponents=4);
    Cryoperm->AddMaterial(Ni, fractionmass=0.775);
    Cryoperm->AddMaterial(Cu, fractionmass=0.045);
    Cryoperm->AddMaterial(Mo, fractionmass=0.025);
    Cryoperm->AddMaterial(Fe, fractionmass=0.155);

    // Define vacuum
    auto *vacuum = man->FindOrBuildMaterial("G4_Galactic");

    // Define the Stainless Steel
    G4double SS_density = 8.03*g/cm3;
    auto *SS = new G4Material("Stainless Steel",density=SS_density,ncomponents=8);
    SS->AddMaterial(N, fractionmass=0.01);
    SS->AddMaterial(Ni, fractionmass=0.1);
    SS->AddMaterial(Cr, fractionmass=0.185);
    SS->AddMaterial(S, fractionmass=0.03);
    SS->AddMaterial(Ph, fractionmass=0.045);
    SS->AddMaterial(Mn, fractionmass=0.02);
    SS->AddMaterial(C, fractionmass=0.08);
    SS->AddMaterial(Fe, fractionmass=0.53);

    // Si3N4
    G4double A;
    G4double Z;
    auto *elSi = new G4Element("Silicon", "Si", Z=14., A=28.0855*g/mole);
    auto *elN = new G4Element("Nitrogen","N", Z=7., A=14.00674*g/mole);
    G4double Si3N4density = 3.44*g/cm3;
    auto *Si3N4 = new G4Material("Si3N4", Si3N4density, 2);
    Si3N4->AddElement (elSi, 3);
    Si3N4->AddElement (elN, 4);

    // -------------------------------------------------------------------------------------------
    // ----------------------------------- Define World volume -----------------------------------
    G4double expHall_x = 3*m;
    G4double expHall_y = 3*m;
    G4double expHall_z = 3*m;

    auto *experimentalHall_box = new G4Box("experimentalHall_box",expHall_x,expHall_y,expHall_z);
    fExperimentalHall_log = new G4LogicalVolume(experimentalHall_box,vacuum,"experimentalHall_log",0,0,0);// vacuum ok, its the mother volume
    fExperimentalHall_phys = new G4PVPlacement(nullptr,G4ThreeVector(), fExperimentalHall_log,"expHall",0,false,0);

    // -------------------------------------------------------------------------------------------
    // ------------------------------------ Define TES array -------------------------------------
    //absorber: 3 µm Bi
    G4double Biabsorberthickness = 3*um;
    G4double Gridthickness = 200*um;
    G4double membrane_thickness = 0.75*um;
    G4double wafer_thickness = 0.25*mm;

    //height of the entire volume to be replicated
    G4double xifu_z = Biabsorberthickness+Gridthickness+membrane_thickness;

    //size of the mother volume to be replicated
    G4double pxl_pitch = 0.814*mm;//V2 res90
    G4double pxl_gap = 11*um;
    G4double pxl_size = pxl_pitch-pxl_gap;
    G4double grid_size = pxl_pitch;

    G4double DetPos_x = 0.0*mm;
    G4double DetPos_y = 0.0*mm;
    G4double DetPos_z = -xifu_z; //So that the absorbers surface is at z=0

    //mother volume of the detector has hexagonal shape
    G4double phiStart = 0;
    G4double phiTotal = 0;
    G4int numSide = 6;
    G4int numZPlanes = 2;
    G4double rInner[2] = {0,0};

    //end of parameters definitions
    auto *element_box = new G4Box("element_box",grid_size*0.5,grid_size*0.5,xifu_z*0.5);
    fElement_log = new G4LogicalVolume(element_box,vacuum,"element_log");

    //definition of the mother volume containing the replicated elements
    G4String pName = "Detector_hex";
    G4double zPlane[2] = {-xifu_z*0.5,xifu_z*0.5};

    G4double rOuter[2] = {8.5*mm,8.5*mm};//apothem

    G4Polyhedra* Detector_hex = new G4Polyhedra(pName,
    phiStart,
    phiTotal,
    numSide,
    numZPlanes,
    zPlane,
    rInner,
    rOuter);

    fDetector_log = new G4LogicalVolume(Detector_hex,vacuum,"Detector_log");//vacuum ok, it's the mother volume
    fDetector_phys = new G4PVPlacement(nullptr, G4ThreeVector(DetPos_x,DetPos_y,DetPos_z+0.0001*um),fDetector_log,"Detector",fExperimentalHall_log, false, 0);

    //Bi absorber 3 um
    auto *Bipxl_box = new G4Box("Bipxl_box",pxl_size*0.5,pxl_size*0.5,Biabsorberthickness*0.5);
    fBipxl_log = new G4LogicalVolume(Bipxl_box,Bi,"Bipxl_log");
    fBipxl_phys = new G4PVPlacement(nullptr,G4ThreeVector(0., 0., xifu_z*0.5-Biabsorberthickness*0.5), fBipxl_log, "Bipxl", fElement_log, false, 0);

    //SiN membrane sustaing the absorbers
    auto *membranepxl_box = new G4Box("membranepxl_box",pxl_size*0.5,pxl_size*0.5,membrane_thickness*0.5);
    fMembranepxl_log = new G4LogicalVolume(membranepxl_box,Si3N4,"membranepxl_log");
    fMembranepxl_phys = new G4PVPlacement(nullptr,G4ThreeVector(0., 0., xifu_z*0.5-Biabsorberthickness-membrane_thickness*0.5), fMembranepxl_log, "membrane", fElement_log, false, 0);

    //step limiter: the membrane is really thin, thinner than the typical interaction step. This forces the particles to perform at least 5 interactions inside
    G4double maxStepmembrane = 0.2*membrane_thickness;
    fStepLimitmembrane = new G4UserLimits(maxStepmembrane);
    fMembranepxl_log->SetUserLimits(fStepLimitmembrane);

    // Grid element
    auto *gridbox = new G4Box("grid_box",grid_size*0.5,grid_size*0.5,Gridthickness*0.5);
    G4double gridbeam_size = 165*um;
    G4double gridhole_size = grid_size-gridbeam_size;
    G4double gridhole_thickness = Gridthickness+1*mm;
    auto *gridhole_box = new G4Box("gridhole_box",gridhole_size*0.5,gridhole_size*0.5,gridhole_thickness*0.5);
    G4SubtractionSolid* gridsubtraction = new G4SubtractionSolid("gridsubtraction",gridbox, gridhole_box);
    fGridpiece_log = new G4LogicalVolume(gridsubtraction,Si,"gridpiece_log");
    fGridpiece_phys = new G4PVPlacement(nullptr,G4ThreeVector(0., 0., xifu_z*0.5-Biabsorberthickness-membrane_thickness-Gridthickness*0.5), fGridpiece_log, "gridpiece", fElement_log, false, 0);

    //parametrization of element_log
    auto *DetParam = new XrayTESdetDetParameterisation(
                317,//V1 and V2 res80 and res90 // to be commented for overlap check
                //2,  // To enable volumes overlap check the number of replicas must be reduced. Keep it commented if not needed
                -grid_size*38,
                grid_size,
                grid_size*0.5,
                grid_size );

    fPhysiDet = new G4PVParameterised(
                  "Det",             // their name
                  fElement_log,       // their logical volume
                  fDetector_log,      // Mother logical volume
                  kXAxis,            // Are placed along this axis
                  317,	             //V1 and V2 res80 and res90
                  //2,               // overlap check. Same as before.
                  DetParam);         // The parametrisation

    // ------------ Define support wafer ----------------------------------------------
    G4double rOutersupp[2] = {39.4*mm,39.4*0.5*mm};
    G4double zwaferPlane[2] = {-0.5*wafer_thickness,0.5*wafer_thickness};

    auto *wafer_hex = new G4Polyhedra("wafer_hex",
    phiStart,
    phiTotal,
    numSide,
    numZPlanes,
    zwaferPlane,
    rOuter,
    rOutersupp);

    fWafer_log = new G4LogicalVolume(wafer_hex,Si,"wafer_log");
    G4double wafer_z = -xifu_z*0.5-Biabsorberthickness-Gridthickness*0.5;
    fWafer_phys = new G4PVPlacement(nullptr, G4ThreeVector(DetPos_x,DetPos_y,wafer_z), fWafer_log,"wafer",fExperimentalHall_log, false, 0);

    //------------------------------------------------------------------------------------------------------
    //----------------------------------   BSC electron detector     ---------------------------------------
    // This fake surface placed just above the absorbers gathers information on the electrons crossing it
    // it can be used to check if the same particle crossed it twice, to identify backscattering events
    //------------------------------------------------------------------------------------------------------
    G4double BSC_z = 0.05*mm;  // half thickness

    G4String pBSCName = "BSCdetector";
    G4double zPlaneBSC[2] = {-BSC_z,BSC_z};
    G4double rOuterBSC[2] = {rOuter[0]+2*mm,rOuter[1]+2*mm};

    auto *BSC_hex = new G4Polyhedra(pBSCName,
      phiStart,
      phiTotal,
      numSide,
      numZPlanes,
      zPlaneBSC,
      rInner,
      rOuterBSC);

    fBSC_log = new G4LogicalVolume(BSC_hex,vacuum,"BSC_log");

    G4double BSCPos_x = 0.0*mm;
    G4double BSCPos_y = 0.0*mm;
    G4double BSCPos_z = wafer_z + wafer_thickness+BSC_z*0.5;

    fBSC_phys = new G4PVPlacement(nullptr, G4ThreeVector(BSCPos_x,BSCPos_y,BSCPos_z),fBSC_log,"BSC",fExperimentalHall_log, false, 0);

    // --------------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------- Define ACD inner array ------------------------------------------------------------
    // Anti Coincidence Detectors are usually thicker than the main ones. It can generate a huge output file if low cuts are used
    // defining a inner part of the ACD and assigning it to a lower precision region, thus allowing only the creation of high energy secondaries
    // can help speeding up the simulation and reducing the output size. In this example only the external 25 um belong to the higher precision region
    // --------------------------------------------------------------------------------------------------------------------------------------------
    G4double ACDthickness = 200*um;
    G4double ACD_distance = 0.75*mm; // distance from the lower edge of the main detector
    G4double ACDarrayPos_z = -1*(ACD_distance+Biabsorberthickness/2+ACDthickness/2);
    G4double ACDpxlGAP = 50*um;

    G4double pxlheight = 8.5*mm;
    G4double pxlwideside = pxlheight*1.191*mm;
    G4double pxlshortside = pxlwideside*0.5;

    G4double inACDthickness = ACDthickness-50*um;

    G4double inpxlheight = pxlheight-50*um;
    G4double inpxlwideside = pxlwideside-50*um;
    G4double inpxlshortside = inpxlwideside*0.5;

    auto *inACDpxl = new G4Trap("inACDpxl",
      inACDthickness,
      inpxlheight,
      inpxlwideside,
      inpxlshortside);

    fInACDpxl_log = new G4LogicalVolume(inACDpxl,Si,"inACDpxl_log");
    G4RotationMatrix rotm = G4RotationMatrix();

    //pxl1
    G4double ACDarrayPos1_x = (pxlshortside*1.5+ACDpxlGAP)*0.5;
    G4double ACDarrayPos1_y = (pxlheight+ACDpxlGAP)*0.5;
    fInACDpxl_phys = new G4PVPlacement(nullptr, G4ThreeVector(ACDarrayPos1_x,ACDarrayPos1_y,ACDarrayPos_z),fInACDpxl_log,"pxl1",fExperimentalHall_log, false, 0);

    //pxl2
    rotm.rotateZ(180*deg);
    G4double ACDarrayPos2_x = -ACDarrayPos1_x;
    G4double ACDarrayPos2_y = -ACDarrayPos1_y;
    G4ThreeVector position2 = G4ThreeVector(ACDarrayPos2_x,ACDarrayPos2_y,ACDarrayPos_z);
    G4Transform3D transform2 = G4Transform3D(rotm,position2);
    fInACDpxl_phys = new G4PVPlacement(transform2,fInACDpxl_log,"pxl2",fExperimentalHall_log, false, 0);

    //pxl3
    rotm.rotateY(180*deg);
    G4double ACDarrayPos3_x = ACDarrayPos1_x;
    G4double ACDarrayPos3_y = -ACDarrayPos1_y;
    G4ThreeVector position3 = G4ThreeVector(ACDarrayPos3_x,ACDarrayPos3_y,ACDarrayPos_z);
    G4Transform3D transform3 = G4Transform3D(rotm,position3);
    fInACDpxl_phys = new G4PVPlacement(transform3,fInACDpxl_log,"pxl3",fExperimentalHall_log, false, 0);

    //pxl4
    rotm.rotateZ(180*deg);
    G4double ACDarrayPos4_x = -ACDarrayPos1_x;
    G4double ACDarrayPos4_y = ACDarrayPos1_y;
    G4ThreeVector position4 = G4ThreeVector(ACDarrayPos4_x,ACDarrayPos4_y,ACDarrayPos_z);
    G4Transform3D transform4 = G4Transform3D(rotm,position4);
    fInACDpxl_phys = new G4PVPlacement(transform4,fInACDpxl_log,"pxl4",fExperimentalHall_log, false, 0);

    // ----------------------------------------------------------------------------------
    // -------------------- Define ACD outer layer array --------------------------------
    auto *ACDpxl = new G4Trap("ACDpxl",
      ACDthickness,
      pxlheight,
      pxlwideside,
      pxlshortside);
    G4SubtractionSolid* ACDsub = new G4SubtractionSolid("ACDsub",ACDpxl, inACDpxl);
    fACDpxl_log = new G4LogicalVolume(ACDsub,Si,"ACDpxl_log");

    //pxl1
    fACDpxl_phys = new G4PVPlacement(nullptr, G4ThreeVector(ACDarrayPos1_x,ACDarrayPos1_y,ACDarrayPos_z),fACDpxl_log,"pxl1a",fExperimentalHall_log, false, 0);

    //pxl2
    fACDpxl_phys = new G4PVPlacement(transform2,fACDpxl_log,"pxl2a",fExperimentalHall_log, false, 0);

    //pxl3
    fACDpxl_phys = new G4PVPlacement(transform3,fACDpxl_log,"pxl3a",fExperimentalHall_log, false, 0);

    //pxl4
    fACDpxl_phys = new G4PVPlacement(transform4,fACDpxl_log,"pxl4a",fExperimentalHall_log, false, 0);

    //------------------------------------------------------------------------------------------------------
    //-----------------------------------------SUPPORTS------------------------------------------------//

    //---------------------   ACD plate     --------------
    G4double ACplate_thickness = 3.5*mm;  // support thickness

    G4String ACplateName = "ACplateName";
    G4double ACplatezPlane[2] = {-0.5*ACplate_thickness,0.5*ACplate_thickness};
    G4double ACplaterInner[2] = {pxlwideside,pxlwideside};//apothem
    G4double ACplaterOuter[2] = {43*mm,43*mm};

    auto *ACplate_hex = new G4Polyhedra(ACplateName,
    phiStart,
    phiTotal,
    numSide,
    numZPlanes,
    ACplatezPlane,
    ACplaterInner,
    ACplaterOuter );

    fACplate_log = new G4LogicalVolume(ACplate_hex,Cu,"ACplate_log");//tbc

    G4double ACplatePos_x = DetPos_x;
    G4double ACplatePos_y = DetPos_y;
    G4double ACplatePos_z = -2.25*mm;

    fACplate_phys = new G4PVPlacement(nullptr, G4ThreeVector(ACplatePos_x,ACplatePos_y,ACplatePos_z),fACplate_log,"ACDplate",fExperimentalHall_log, false, 0);

    //------------------------------------------------------------------------------------------------------
    //----------------------------------------------   cage   ----------------------------------------------
    // ALL THE FOLLOWING SUPPORTING STRUCTURES DEFINE A COMPLEX SHAPE THAT WAS SIMPLIFIED AND ADAPTED FROM AN ENGINEERING CAD MODEL
    // ALL THE NUMBERS RELATIVE TO SIZES AND POSITIONS HAVE BEEN EXTRACTED FROM THE CAD FILE, THAT IS WHY THE NUMBERS ARE SO ODD
    // THIS IS MEANT JUST TO GIVE AN EXAMPLE OF THE SUPPORTING STRUCTURES PRESENT IN AN ACTUAL CRYOSTAT

    //------------------supporting columns (full)
    G4double hightOfcagecolumn = 60*mm;

    auto *cagecolumn = new G4Tubs("cagecolumn", 0*mm, 3.99*mm, hightOfcagecolumn*0.5, 0*deg, 360*deg);
    G4LogicalVolume* cagecolumn_log = new G4LogicalVolume(cagecolumn,Cu,"cagecolumn_log",0,0,0);

    fCagecolumn_phys = new G4PVPlacement(
      0,							                          // no rotation
      G4ThreeVector(45.75*mm,0*mm,-36.5027*mm),	// translation position
      cagecolumn_log,					                  // its logical volume
      "cagecolumn1",                             // its name
      fExperimentalHall_log,					            // its mother (logical) volume
      false,							                      // no boolean operations
      0);

    fCagecolumn2_phys = new G4PVPlacement(
      0,
      G4ThreeVector(-45.75*mm,0*mm,-36.5027*mm),
      cagecolumn_log,
      "cagecolumn2",
      fExperimentalHall_log,
      false,
      0);

    fCagecolumn3_phys = new G4PVPlacement(
      0,
      G4ThreeVector(22.875*mm,39.621*mm,-36.5027*mm),
      cagecolumn_log,
      "cagecolumn3",
      fExperimentalHall_log,
      false,
      0);

    fCagecolumn4_phys = new G4PVPlacement(
      0,
      G4ThreeVector(22.875*mm,-39.621*mm,-36.5027*mm),
      cagecolumn_log,
      "cagecolumn4",
      fExperimentalHall_log,
      false,
      0);

    fCagecolumn5_phys = new G4PVPlacement(
      0,
      G4ThreeVector(-22.875*mm,-39.621*mm,-36.5027*mm),
      cagecolumn_log,
      "cagecolumn5",
      fExperimentalHall_log,
      false,
      0);

    fCagecolumn6_phys = new G4PVPlacement(
      0,
      G4ThreeVector(-22.875*mm,39.621*mm,-36.5027*mm),
      cagecolumn_log,
      "cagecolumn6",
      fExperimentalHall_log,
      false,
      0);

    //----------------------------------supporting columns (empty)
    auto *emptycagecolumn = new G4Tubs("emptycagecolumn", 1.25*mm, 3*mm, hightOfcagecolumn*0.5, 0*deg, 360*deg);
    G4LogicalVolume* emptycagecolumn_log = new G4LogicalVolume(emptycagecolumn,Cu,"emptycagecolumn_log",0,0,0);

    fEmptycagecolumn_phys = new G4PVPlacement(
      0,
      G4ThreeVector(22.875*mm,0*mm,-36.5027*mm),
      emptycagecolumn_log,
      "emptycagecolumn1",
      fExperimentalHall_log,
      false,
      0);

    fEmptycagecolumn2_phys = new G4PVPlacement(
      0,
      G4ThreeVector(-11.4375*mm,-19.8103311*mm,-36.5027*mm),
      emptycagecolumn_log,
      "emptycagecolumn2",
      fExperimentalHall_log,
      false,
      0);

    fEmptycagecolumn3_phys = new G4PVPlacement(
      0,
      G4ThreeVector(-11.4375*mm,19.810331*mm,-36.5027*mm),
      emptycagecolumn_log,
      "emptycagecolumn3",
      fExperimentalHall_log,
      false,
      0);

    //-----------------------walls
    G4double cagewall_x = 1*mm;
    G4double cagewall_y = 32.62*mm;
    G4double cagewall_z = 65*mm;

    auto *cagewall_box = new G4Box("_box",cagewall_x*0.5,cagewall_y*0.5,cagewall_z*0.5);

    fCagewall_log = new G4LogicalVolume(cagewall_box,Cu,"cagewall_log",0,0,0);  //tbc
    fCagewall_phys = new G4PVPlacement(
      0,
      G4ThreeVector(22.875*mm,19.3103311*mm,-36.5027*mm),
      fCagewall_log,
      "cagewall1",
      fExperimentalHall_log,
      false,
      0);

    ///and the other 5 pieces
    fCagewall2_phys = new G4PVPlacement(
      0,
      G4ThreeVector(22.875*mm,-19.3103311*mm,-36.5027*mm),
      fCagewall_log,
      "cagewall2",
      fExperimentalHall_log,
      false,
      0);

    G4RotationMatrix wallrotm = G4RotationMatrix();

    wallrotm.rotateZ(-60*deg);
    G4ThreeVector wallposition3 = G4ThreeVector(5.2857373*mm,29.4654967*mm,-36.5027*mm);
    G4Transform3D walltransform3 = G4Transform3D(wallrotm,wallposition3);
    fCagewall3_phys = new G4PVPlacement(walltransform3,fCagewall_log,"cagewall3",fExperimentalHall_log, true, 0);

    G4ThreeVector wallposition4 = G4ThreeVector(-28.1607373*mm,10.1551656*mm,-36.5027*mm);
    G4Transform3D walltransform4 = G4Transform3D(wallrotm,wallposition4);
    fCagewall4_phys = new G4PVPlacement(walltransform4,fCagewall_log,"cagewall4",fExperimentalHall_log, true, 0);

    wallrotm.rotateZ(120*deg);
    G4ThreeVector wallposition5 = G4ThreeVector(-28.1607373*mm,-10.1551656*mm,-36.5027*mm);
    G4Transform3D walltransform5 = G4Transform3D(wallrotm,wallposition5);
    fCagewall5_phys = new G4PVPlacement(walltransform5,fCagewall_log,"cagewall5",fExperimentalHall_log, true, 0);

    G4ThreeVector wallposition6 = G4ThreeVector(5.2857373*mm,-29.4654967*mm,-36.5027*mm);
    G4Transform3D walltransform6 = G4Transform3D(wallrotm,wallposition6);
    fCagewall6_phys = new G4PVPlacement(walltransform6,fCagewall_log,"cagewall6",fExperimentalHall_log, true, 0);

    //-----------------------miniwalls
    G4double minicagewall_x = 15.875*mm;
    G4double minicagewall_y = 1*mm;

    auto *minicagewall_box = new G4Box("minicagewall_box",minicagewall_x*0.5,minicagewall_y*0.5,cagewall_z*0.5);
    fMinicagewall_log = new G4LogicalVolume(minicagewall_box,Cu,"minicagewall_log",0,0,0);//tbc
    fMinicagewall_phys = new G4PVPlacement(
      0,
      G4ThreeVector(33.8125*mm,0*mm,-36.5027*mm),
      fMinicagewall_log,
      "minicagewall1",
      fExperimentalHall_log,
      false,
      0);

    G4RotationMatrix wallrotm2 = G4RotationMatrix();

    wallrotm2.rotateZ(-60*deg);
    G4ThreeVector miniwallposition2 = G4ThreeVector(-16.90625*mm,29.282484*mm,-36.5027*mm);
    G4Transform3D miniwalltransform2 = G4Transform3D(wallrotm2,miniwallposition2);
    fMinicagewall2_phys = new G4PVPlacement(miniwalltransform2,fMinicagewall_log,"minicagewall2",fExperimentalHall_log, true, 0);

    wallrotm2.rotateZ(-60*deg);
    G4ThreeVector miniwallposition3 = G4ThreeVector(-16.90625*mm,-29.282484*mm,-36.5027*mm);
    G4Transform3D miniwalltransform3 = G4Transform3D(wallrotm2,miniwallposition3);
    fMinicagewall3_phys = new G4PVPlacement(miniwalltransform3,fMinicagewall_log,"minicagewall3",fExperimentalHall_log, true, 0);

    //--------------extwalls
    G4double extcagewall_x = 1.85*mm;
    G4double extcagewall_y = 37.75*mm;
    G4double extcagewall_z = 65*mm;

    auto *extcagewall_box = new G4Box("extcagewall_box",extcagewall_x*0.5,extcagewall_y*0.5,extcagewall_z*0.5);
    fExtcagewall_log = new G4LogicalVolume(extcagewall_box,Cu,"extcagewall_log",0,0,0);
    G4RotationMatrix wallrotm3 = G4RotationMatrix();
    wallrotm3.rotateZ(-30*deg);
    G4ThreeVector extwallposition = G4ThreeVector(-35.0724512*mm,20.2490891*mm,-36.5027*mm);
    G4Transform3D extwalltransform = G4Transform3D(wallrotm3,extwallposition);
    fExtcagewall_phys = new G4PVPlacement(extwalltransform,fExtcagewall_log,"extcagewall1",fExperimentalHall_log, true, 0);

    G4ThreeVector extwallposition2 = G4ThreeVector(35.0724512*mm,-20.2490891*mm,-36.5027*mm);
    G4Transform3D extwalltransform2 = G4Transform3D(wallrotm3,extwallposition2);
    fExtcagewall2_phys = new G4PVPlacement(extwalltransform2,fExtcagewall_log,"extcagewall2",fExperimentalHall_log, true, 0);

    //--------------extboards
    G4double extboard_x = 0.5*mm;
    G4double extboard_y = 36.5*mm;
    G4double extboard_z = 64.5*mm;

    auto *extboard_box = new G4Box("extboard_box",extboard_x*0.5,extboard_y*0.5,extboard_z*0.5);
    fExtboard_log = new G4LogicalVolume(extboard_box,Si,"extboard_log",0,0,0);
    G4ThreeVector extboardposition = G4ThreeVector(37.932*mm,-21.9*mm,-36.253*mm);
    G4Transform3D extboardtransform = G4Transform3D(wallrotm3,extboardposition);
    fExtboard_phys = new G4PVPlacement(extboardtransform,fExtboard_log,"extboard1",fExperimentalHall_log, true, 0);

    G4ThreeVector extboardposition2 = G4ThreeVector(-37.932*mm,21.9*mm,-36.253*mm);
    G4Transform3D extboardtransform2 = G4Transform3D(wallrotm3,extboardposition2);
    fExtboard2_phys = new G4PVPlacement(extboardtransform2,fExtboard_log,"extboard2",fExperimentalHall_log, true, 0);

    //-------------------------------------------------------------------------------------------------
    // --- END OF THE SUPPORTING STRUCTURES. HERE: STAGES OF THE CRYOSTAT -----------------------------

    //-------------------------------------------------------------------------------------------------
    //----------------------------------------   First thermal shield  -----------------------------------------

    auto *FirstShield_botplate = new G4Tubs("Tube_FirstShield_botplate", 0.*mm, 80*mm, 1*mm*0.5, 0*deg, 360*deg);
    G4LogicalVolume* FirstShield_botplate_log = new G4LogicalVolume(FirstShield_botplate,Cryoperm,"FirstShield_botplate_log",0,0,0);
    fFirstShield_botplate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-101*mm),
      FirstShield_botplate_log,
      "FirstShield_botplate",
      fExperimentalHall_log,
      false,
      0);

    //------------------------
    auto *FirstShield_side = new G4Tubs("Tube_FirstShield_side", 79*mm, 80*mm, 145*mm*0.5, 0*deg, 360*deg);
    G4LogicalVolume* FirstShield_side_log = new G4LogicalVolume(FirstShield_side,Cryoperm,"FirstShield_side_log",0,0,0);
    fFirstShield_side_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-28*mm),
      FirstShield_side_log,
      "FirstShield_side",
      fExperimentalHall_log,
      false,
      0);

    //---------------------------
    auto *FirstShield_topplate = new G4Tubs("Tube_FirstShield_topplate", 31.85*mm, 80*mm, 1*mm*0.5, 0*deg, 360*deg);
    G4LogicalVolume* FirstShield_topplate_log = new G4LogicalVolume(FirstShield_topplate,Cryoperm,"FirstShield_topplate_log",0,0,0);
    fFirstShield_topplate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,45*mm),
      FirstShield_topplate_log,
      "FirstShield_topplate",
      fExperimentalHall_log,
      false,
      0);

    //------------------------------
    auto *FirstShield_topcyl = new G4Tubs("Tube_FirstShield_topcyl", 31.85*mm, 32.85*mm, 74*mm*0.5, 0*deg, 360*deg);
    G4LogicalVolume* FirstShield_topcyl_log = new G4LogicalVolume(FirstShield_topcyl,Cryoperm,"FirstShield_topcyl_log",0,0,0);
    fFirstShield_topcyl_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,82.5*mm),
      FirstShield_topcyl_log,
      "FirstShield_topcyl",
      fExperimentalHall_log,
      false,
      0);

    //-------------------------------------------------------------------------------------------------
    //------------------------------------------   Ring  ----------------------------------------------

    // This solid has a complex shape, filled with holes and supporting beams. Since it is not in direct sight of
    // the detector, it can be simplified as a ring-like shape.

    // The material has been defined as a custom Al with density equal to the average density in the original solid
    G4double hightOfRing = 85*mm;
    auto *Ring = new G4Tubs("Ring_tubs", 80*mm, 97.1*mm, hightOfRing*0.5, 0*deg, 360*deg);

    // Three boxes are used for a boolean subtraction, cutting 3 sections of the ring
    G4double cuts_x = 38.781347139*mm;
    G4double cuts_y = 116.620103593*mm;
    G4double cuts_z = 85.1*mm;

    auto *cuts_box = new G4Box("cuts_box",cuts_x*0.5,cuts_y*0.5,cuts_z*0.5);

    G4ThreeVector cutspos = G4ThreeVector(-105.055643631*mm,0,0);
    G4SubtractionSolid* Ringsub = new G4SubtractionSolid("Ringsub",Ring,cuts_box,0,cutspos);

    auto *cutsrotm = new G4RotationMatrix();
    cutsrotm->rotateZ(60*deg);
    G4ThreeVector cutspos2 = G4ThreeVector(52.527821885*mm,-90.980856208*mm,0);
    G4SubtractionSolid* Ringsub2 = new G4SubtractionSolid("Ringsub2",Ringsub,cuts_box,cutsrotm,cutspos2);
    cutsrotm->rotateZ(60*deg);
    G4ThreeVector cutspos3 = G4ThreeVector(52.527821885*mm,90.980856208*mm,0);
    auto *Ringsub3 = new G4SubtractionSolid("Ringsub3",Ringsub2,cuts_box,cutsrotm,cutspos3);
    auto *Ring_log = new G4LogicalVolume(Ringsub3,RingAl,"Ring_log",0,0,0);//tbc
    fRing_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,-22*mm),
      Ring_log,
      "Ring",
      fExperimentalHall_log,
      false,
      0);

    //-------------------------------------------------------------------------------------------------
    //----------------------------------------   Second thermal shield  -----------------------------------------
    // This is another irregular and complex shape derived from a CAD engineering model that was simplified for Geant4 implementation

    // Define the overall shape and bottom plate defining the base points to extrude
    std::vector<G4TwoVector> newShape(6);
    newShape[0] = G4TwoVector(58.936*mm, 154.115*mm);
    newShape[1] = G4TwoVector(104*mm, 128.098*mm);
    newShape[2] = G4TwoVector(104*mm, -128.098*mm);
    newShape[3] = G4TwoVector(58.936*mm, -154.115*mm);
    newShape[4] = G4TwoVector(-162.936*mm, -26.018*mm);
    newShape[5] = G4TwoVector(-162.936*mm, 26.018*mm);

    // ------------------------------------------- the lower plate  -------------------------------------------
    // Extrusion of the solid defined by the 6 points above
    auto *MyShape = new G4ExtrudedSolid("MyShape", newShape, 2.65*0.5*mm, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);

    // Inserting a hole in the plate
    auto *MyHole = new G4Tubs("MyHole", 0, 19*mm, 3*mm*0.5, 0*deg, 360*deg);
    G4ThreeVector holepos = G4ThreeVector(0,0,0);
    auto *SecondShieldbotPlate = new G4SubtractionSolid("SecondShieldbotPlate_sub",MyShape,MyHole,0,holepos);
    auto *SecondShieldbotPlate_log = new G4LogicalVolume(SecondShieldbotPlate,Al,"SecondShieldbotPlate_log",0,0,0);
    fSecondShieldbotPlate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-110.178*mm),
      SecondShieldbotPlate_log,
      "SecondShieldbotPlate",
      fExperimentalHall_log,
      false,
      0);

    //---------------------------------------------------------------------------------------------------
    //------------------------------------------ the central part ---------------------------------------

    // Definitino of the central block
    auto *MyShape2 = new G4ExtrudedSolid("MyShape2", newShape, 215.35*0.5*mm, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);

    // The central part needs to be empty inside, so we define 6 new points for the subtraction of the central area
    std::vector<G4TwoVector> newHole(6);
    newHole[0] = G4TwoVector(58.936*mm, 149.762*mm);
    newHole[1] = G4TwoVector(100.23*mm, 125.921*mm);
    newHole[2] = G4TwoVector(100.23*mm, -125.921*mm);
    newHole[3] = G4TwoVector(58.936*mm, -149.762*mm);
    newHole[4] = G4TwoVector(-159.166*mm, -23.841*mm);
    newHole[5] = G4TwoVector(-159.166*mm, 23.841*mm);

    // Extrusion of the subtraction solid
    auto *MyInnerShape = new G4ExtrudedSolid("MyInnerShape", newHole, 216*0.5*mm, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);

    // creation of the central part: subtract the  subtraction solid from the central block
    auto *SecondShieldside = new G4SubtractionSolid("SecondShieldside_sub",MyShape2,MyInnerShape,0,holepos);
    auto *SecondShieldside_log = new G4LogicalVolume(SecondShieldside,Al,"SecondShieldside_log",0,0,0);
    fSecondShieldside_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-1.178*mm),
      SecondShieldside_log,
      "SecondShieldside",
      fExperimentalHall_log,
      false,
      0);

    //-------------------------------------- top plate cover ---------------------------------------
    G4ExtrudedSolid *MyShape3 = new G4ExtrudedSolid("MyShape3", newShape, 10.87*0.5*mm, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);
    auto *MyHole2 = new G4Tubs("MyHole2", 0, 36.85*mm, 11*mm*0.5, 0*deg, 360*deg);
    auto *SecondShieldtopPlate = new G4SubtractionSolid("SecondShieldtopPlate_sub",MyShape3,MyHole2,0,holepos);
    auto *SecondShieldtopPlate_log = new G4LogicalVolume(SecondShieldtopPlate,Al,"SecondShieldtopPlate_log",0,0,0);
    fSecondShieldtopPlate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,111.932*mm),
      SecondShieldtopPlate_log,
      "SecondShieldtopPlate",
      fExperimentalHall_log,
      false,
      0);

    //--------------------------------------
    auto *SecondShieldshield_topcyl = new G4Tubs("SecondShieldshield", 36.85*mm, 39*mm, 12*mm*0.5, 0*deg, 360*deg);
    auto *SecondShieldshield_topcyl_log = new G4LogicalVolume(SecondShieldshield_topcyl,Al,"SecondShieldshield_topcyl_log",0,0,0);
    fSecondShieldshield_topcyl_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,123.367*mm),
      SecondShieldshield_topcyl_log,
      "SecondShieldshield_topcyl",
      fExperimentalHall_log,
      false,
      0);

    //------------------------------------------------------------------------------------------------------
    //-------------------------------------    EXTERNAL ZONE     -------------------------------------------

    //------------------------------------------------------------------------------------------------------
    //---------------------------------------- Third thermal shield ---------------------------------------------------
    G4double TopConeDOWNradius = 174.6*mm;
    G4double TopConeUPradius = 93*mm;
    G4double TopConeheight = 40*mm;
    G4double SideHeight = 280*mm;
    G4double ConeDOWNradius = 80*mm;
    G4double Conethick = 0.4*mm;
    G4double ConeUPradius = 174.6*mm;
    G4double Coneheight = 40*mm;

    //---------------
    auto *ThirdShieldside = new G4Tubs("ThirdShieldside_tubs", TopConeDOWNradius, TopConeDOWNradius+Conethick, SideHeight*0.5, 0*deg, 360*deg);
    G4LogicalVolume* ThirdShieldside_log = new G4LogicalVolume(ThirdShieldside,Al,"ThirdShieldside_log",0,0,0);
    fThirdShieldside_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,0),
      ThirdShieldside_log,
      "ThirdShieldside",
      fExperimentalHall_log,
      false,
      0);

    //-------------------
    auto *ThirdShieldtopCone = new G4Cons("ThirdShieldtopCone_cons", TopConeDOWNradius,Conethick+TopConeDOWNradius,TopConeUPradius,TopConeUPradius+Conethick, TopConeheight*0.5, 0*deg, 360*deg);
    auto *ThirdShieldtopCone_log = new G4LogicalVolume(ThirdShieldtopCone,Al,"ThirdShieldtopCone_log",0,0,0);
    fThirdShieldtopCone_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,(SideHeight+TopConeheight)*0.5),
      ThirdShieldtopCone_log,
      "ThirdShieldtopCone",
      fExperimentalHall_log,
      false,
      0);

    //-------------------
    auto *ThirdShieldbotCone = new G4Cons("ThirdShieldbotCone_cons", ConeDOWNradius,Conethick+ConeDOWNradius,ConeUPradius,ConeUPradius+Conethick, Coneheight*0.5, 0*deg, 360*deg);
    auto *ThirdShieldbotCone_log = new G4LogicalVolume(ThirdShieldbotCone,Al,"ThirdShieldbotCone_log",0,0,0);
    fThirdShieldbotCone_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,(-SideHeight-Coneheight)*0.5),
      ThirdShieldbotCone_log,
      "ThirdShieldbotCone",
      fExperimentalHall_log,
      false,
      0);

    //-------------------
    auto *ThirdShieldbotPlate = new G4Tubs("ThirdShieldbotPlate_tubs", 0*mm, Conethick+ConeDOWNradius,Conethick*0.5, 0*deg, 360*deg);
    auto *ThirdShieldbotPlate_log = new G4LogicalVolume(ThirdShieldbotPlate,Al,"ThirdShieldbotPlate_log",0,0,0);
    fThirdShieldbotPlate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,(-SideHeight-Conethick)*0.5-Coneheight),
      ThirdShieldbotPlate_log,
      "ThirdShieldbotPlate",
      fExperimentalHall_log,
      false,
      0);

    //----------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------- External layer of the Cryostat --------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------
    G4double sideradius = 200*mm;
    G4double midsidethick = 8.8*mm;
    G4double midsideheigth = 120*mm;
    auto *CryostatmidSide = new G4Tubs("CryostatmidSide_tubs", sideradius, sideradius+midsidethick, midsideheigth*0.5, 0*deg, 360*deg);
    auto *CryostatmidSide_log = new G4LogicalVolume(CryostatmidSide,Al,"CryostatmidSide_log",0,0,0);
    fCryostatmidSide_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,0),
      CryostatmidSide_log,
      "CryostatmidSide",
      fExperimentalHall_log,
      false,
      0);

    //---------------
    G4double topsidethick = 6.5*mm;
    G4double topsideheigth = 85*mm;
    auto *CryostattopSide = new G4Tubs("CryostattopSide_tubs", sideradius, sideradius+topsidethick, topsideheigth*0.5, 0*deg, 360*deg);
    auto *CryostattopSide_log = new G4LogicalVolume(CryostattopSide,Al,"CryostattopSide_log",0,0,0);
    fCryostattopSide_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,midsideheigth*0.5+topsideheigth*0.5),
      CryostattopSide_log,
      "CryostattopSide",
      fExperimentalHall_log,
      false,
      0);

    //-------------------
    ConeUPradius = 97*mm;
    Coneheight = 75*mm;
    auto *CryostattopCone = new G4Cons("CryostattopCone_cons", sideradius,topsidethick+sideradius,ConeUPradius,ConeUPradius+topsidethick, Coneheight*0.5, 0*deg, 360*deg);
    auto *CryostattopCone_log = new G4LogicalVolume(CryostattopCone,Al,"CryostattopCone_log",0,0,0);
    fCryostattopCone_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,midsideheigth*0.5+topsideheigth+Coneheight*0.5),	// translation position
      CryostattopCone_log,
      "CryostattopCone",
      fExperimentalHall_log,
      false,
      0);

    //---------------
    G4double botsideheigth = 85*mm;
    auto *CryostatbotSide = new G4Tubs("CryostatbotSide_tubs", sideradius, sideradius+topsidethick, botsideheigth*0.5, 0*deg, 360*deg);
    auto *CryostatbotSide_log = new G4LogicalVolume(CryostatbotSide,Al,"CryostatbotSide_log",0,0,0);
    fCryostatbotSide_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-midsideheigth*0.5-botsideheigth*0.5),
      CryostatbotSide_log,
      "CryostatbotSide",
      fExperimentalHall_log,
      false,
      0);

    //-------------------
    ConeDOWNradius = 120*mm;
    auto *CryostatbotCone = new G4Cons("CryostatbotCone_cons", ConeDOWNradius,topsidethick+ConeDOWNradius,sideradius,sideradius+topsidethick, Coneheight*0.5, 0*deg, 360*deg);
    auto *CryostatbotCone_log = new G4LogicalVolume(CryostatbotCone,Al,"CryostatbotCone_log",0,0,0);
    fCryostatbotCone_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,-midsideheigth*0.5-botsideheigth-Coneheight*0.5),	// translation position
      CryostatbotCone_log,
      "CryostatbotCone",
      fExperimentalHall_log,
      false,
      0);

    //---------------------
    auto *CryostatbotPlate = new G4Tubs("CryostatbotPlate_tubs", 0*mm, ConeDOWNradius, midsidethick*0.5, 0*deg, 360*deg);
    auto *CryostatbotPlate_log = new G4LogicalVolume(CryostatbotPlate,Al,"CryostatbotPlate_log",0,0,0);
    fCryostatbotPlate_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0,0,-midsideheigth*0.5-botsideheigth-Coneheight+midsidethick*0.5),	// translation position
      CryostatbotPlate_log,
      "CryostatbotPlate",
      fExperimentalHall_log,
      false,
      0);

    //------------------------------------------------------------------------------------------------------
    //-------------------------------------  radiation filter   --------------------------------------------
    //------------------------------------------------------------------------------------------------------
    G4double halfhightOfAlFilter = 20*0.5*nm;
    G4double mesh_thickness = 80*um;
    G4double filtercarrierthickness = 5*mm;
    G4double outerRadiusOfFilter = 128*0.5*mm;
    G4double ZposFilter = 180*mm;

    auto *AlFilter = new G4Tubs("Filter", 0, outerRadiusOfFilter, halfhightOfAlFilter, 0.*deg, 360.*deg);
    auto *AlFilter_log = new G4LogicalVolume(AlFilter,Al,"AlFilter_log",0,0,0);
    fAlFilter_phys = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,ZposFilter), AlFilter_log, "AlFilter" ,fExperimentalHall_log,false,0);

    //------step limiter filter
    AlFilter_log->SetUserLimits(fStepLimitAlfilta);

    //----- mesh: Nb
    auto *Mesh = new G4Tubs("meshcage", 0, outerRadiusOfFilter+3.3*mm, mesh_thickness*0.5 , 0.*deg, 360.*deg);
    fMesh_log = new G4LogicalVolume(Mesh,vacuum,"Mesh_log");   //vacuum, it's the mother volume and not the grid
    G4double meshPos_z = ZposFilter-halfhightOfAlFilter-mesh_thickness*0.5;
    fMesh_phys = new G4PVPlacement(nullptr, G4ThreeVector(DetPos_x,DetPos_y,meshPos_z),fMesh_log,"meshphys",fExperimentalHall_log, false, 0);

    //////// mesh pixels (grid support)
    G4double half_x1 = 4827/2.*um;
    G4double half_x2 = 4827/2.*um;
    G4double half_x3 = 4880/2.*um;
    G4double half_x4 = 4880/2.*um;
    G4double half_y1 = 80/2.*um;
    G4double half_y2 = 80/2.*um;
    G4double dZ = 26.50/2.*um;

    G4double L1 = half_x1;
    G4double L3 = dZ;

    // create the trapezoids for the element of the grid
    auto *trapezoid_A = new G4Trap("Trapezoid_A",dZ,0,0,half_y1,half_x1,half_x2,0,half_y2,half_x3,half_x4,0);
    auto *trapezoid_B = new G4Trap("Trapezoid_B",dZ,0,0,half_y1,half_x1,half_x2,0,half_y2,half_x3,half_x4,0);
    auto *trapezoid_C = new G4Trap("Trapezoid_C",dZ,0,0,half_y1,half_x1,half_x2,0,half_y2,half_x3,half_x4,0);
    auto *trapezoid_D = new G4Trap("Trapezoid_D",dZ,0,0,half_y1,half_x1,half_x2,0,half_y2,half_x3,half_x4,0);

    // create their rotation matrices and placement vectors
    auto *rotmat_A = new G4RotationMatrix();
    auto *rotmat_B = new G4RotationMatrix();
    auto *rotmat_C = new G4RotationMatrix();
    auto *rotmat_D = new G4RotationMatrix();
    rotmat_A->rotateX(90.*deg);
    rotmat_B->rotateX(90.*deg);
    rotmat_B->rotateY(-90.*deg);
    rotmat_C->rotateX(270.*deg);
    rotmat_D->rotateX(90.*deg);
    rotmat_D->rotateY(90.*deg);
    G4ThreeVector pos_trap_A = G4ThreeVector(0,L1+L3,meshPos_z);
    G4ThreeVector pos_trap_B = G4ThreeVector(L1+L3,0,meshPos_z);
    G4ThreeVector pos_trap_C = G4ThreeVector(0,-L1-L3,meshPos_z);
    G4ThreeVector pos_trap_D = G4ThreeVector(-L1-L3,0,meshPos_z);

    //Test: create their logical and physical volumes
    fTrapezoid_A_log = new G4LogicalVolume(trapezoid_A,SS,"trapezoid_A_log");
    fTrapezoid_B_log = new G4LogicalVolume(trapezoid_B,SS,"trapezoid_B_log");
    fTrapezoid_C_log = new G4LogicalVolume(trapezoid_C,SS,"trapezoid_C_log");
    fTrapezoid_D_log = new G4LogicalVolume(trapezoid_D,SS,"trapezoid_D_log");

    G4ThreeVector pos_munion_base = G4ThreeVector(0,0,0);
    auto *rotmat_final = new G4RotationMatrix();
    rotmat_final->rotateX(0.*deg);
    G4Transform3D transf_base = G4Transform3D(*rotmat_final, pos_munion_base);

    G4int index = 0;
    G4double coordx = 0;
    G4double coordy = 0;
    const G4int num = 13;
    G4double distance = 0;

    G4ThreeVector pos_base = G4ThreeVector(0,0,0);
    G4ThreeVector posA = G4ThreeVector(0,0,0); // the z coord is 0 because it will refer to Mesh_log
    G4ThreeVector posB = G4ThreeVector(0,0,0);
    G4ThreeVector posC = G4ThreeVector(0,0,0);
    G4ThreeVector posD = G4ThreeVector(0,0,0);

    std::vector<G4VPhysicalVolume*> phys_array;
    for (int i=-num; i<=num; i++)
    {
      for (int j=-num; j<=num; j++)
      {
        std::string sA = "trapezoidA_";
        std::string sB = "trapezoidB_";
        std::string sC = "trapezoidC_";
        std::string sD = "trapezoidD_";

        std::string index_str = std::to_string(index);
        sA += index_str;
        sB += index_str;
        sC += index_str;
        sD += index_str;

        coordx = i*half_x3*2;
        coordy = j*half_x3*2;

        distance = std::sqrt(coordx*coordx + coordy*coordy);

        if (distance < outerRadiusOfFilter)
        {
          pos_base = G4ThreeVector(coordx,coordy,meshPos_z);
          posA = G4ThreeVector(coordx,coordy+L1+L3,0); // the z coord is 0 because it will refer to fMesh_log
          posB = G4ThreeVector(coordx+L1+L3,coordy,0);
          posC = G4ThreeVector(coordx,coordy-L1-L3,0);
          posD = G4ThreeVector(coordx-L1-L3,coordy,0);

          fTrapezoid_A_phys = new G4PVPlacement(rotmat_A,posA,fTrapezoid_A_log,sA,fMesh_log, false, 0);
          fTrapezoid_B_phys = new G4PVPlacement(rotmat_B,posB,fTrapezoid_B_log,sB,fMesh_log, false, 0);
          fTrapezoid_C_phys = new G4PVPlacement(rotmat_C,posC,fTrapezoid_C_log,sC,fMesh_log, false, 0);
          fTrapezoid_D_phys = new G4PVPlacement(rotmat_D,posD,fTrapezoid_D_log,sD,fMesh_log, false, 0);

          phys_array.push_back(fTrapezoid_A_phys);
          phys_array.push_back(fTrapezoid_B_phys);
          phys_array.push_back(fTrapezoid_C_phys);
          phys_array.push_back(fTrapezoid_D_phys);
          index += 1;
        }
      }
    }

    auto *filtercarrier = new G4Tubs("filtercarrier_tubs", outerRadiusOfFilter,outerRadiusOfFilter+3.85*mm ,filtercarrierthickness*0.5,0.*deg,360.*deg);
    auto *filtercarrier_log = new G4LogicalVolume(filtercarrier,Al,"filtercarrier_log",0,0,0);
    fFiltercarrier_phys = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,meshPos_z+mesh_thickness*0.5+filtercarrierthickness*0.5),filtercarrier_log,"filtercarrier",fExperimentalHall_log,false,0);

    //------------------------------------------------------------------------------------------------------
    //----------------------------------------   aperture cylinder   ---------------------------------------
    //------------------------------------------------------------------------------------------------------
    //300K piece
    //G4double ACinnerDOWNradius = 106*mm;
    G4double ACinnerUPradius = 67.8644*mm;
    G4double ACthick = 1*mm;
    G4double ACheight = midsideheigth*0.5+topsideheigth+Coneheight-ZposFilter-filtercarrierthickness;

    //ZposFilter
    auto *AC = new G4Cons("AC_cons",ACinnerUPradius-ACthick,ACinnerUPradius,ConeUPradius-ACthick,ConeUPradius,ACheight*0.5,0*deg,360*deg);
    auto *AC_log = new G4LogicalVolume(AC,Ti,"AC_log",0,0,0);
    fAC_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,midsideheigth*0.5+topsideheigth+Coneheight-ACheight*0.5),
      AC_log,
      "AC",
      fExperimentalHall_log,
      false,
      0);

    //-----inner coating
    auto *ACIC = new G4Cons("ACIC_cons",ACinnerUPradius-ACthick-5*um,ACinnerUPradius-ACthick,ConeUPradius-ACthick-5*um,ConeUPradius-ACthick,ACheight*0.5,0*deg,360*deg);
    auto *ACIC_log = new G4LogicalVolume(ACIC,Au,"ACIC_log",0,0,0);
    fACIC_phys = new G4PVPlacement(
      0,
      G4ThreeVector(0.,0.,midsideheigth*0.5+topsideheigth+Coneheight-ACheight*0.5),
      ACIC_log,
      "ACIC",
      fExperimentalHall_log,
      false,
      0);

    //------------------------------------------------------------------------------------------------------
    //------------------------------    The Al sphere representing the spacecraft     ----------------------
    //------------------------------------------------------------------------------------------------------
    auto *sphere = new G4Sphere("sphere_sphere",26*cm,27*cm,0.,twopi,0.,pi);
    fSphere_log = new G4LogicalVolume(sphere,Al,"sphere_log");
    fSphere_phys = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,0*mm),fSphere_log,"sphere",fExperimentalHall_log,false,0);

    //------------------------------------------------------------------------------------------------------
    //---------------------    regions assignment: closer to the detector, lower the Ecut     --------------
    //------------------------------------------------------------------------------------------------------
    // everything directly seen by the detector
    InnerRegion->AddRootLogicalVolume(fElement_log);
    InnerRegion->AddRootLogicalVolume(fDetector_log);
    InnerRegion->AddRootLogicalVolume(fBipxl_log);
    InnerRegion->AddRootLogicalVolume(fMembranepxl_log);
    InnerRegion->AddRootLogicalVolume(fGridpiece_log);
    InnerRegion->AddRootLogicalVolume(fWafer_log);
    InnerRegion->AddRootLogicalVolume(fACDpxl_log);
    InnerRegion->AddRootLogicalVolume(AlFilter_log);  // it was in the intermediate region
    InnerRegion->AddRootLogicalVolume(ACIC_log);

    fWorld_phys = fExperimentalHall_phys;
  }
  else
  {
    G4cout << "Build geometry from GDML file: " << fReadFile << G4endl;
    G4VPhysicalVolume* fWorldPhysVol;
    fParser.Read(fReadFile);
    fWorldPhysVol = fParser.GetWorldVolume();
    fWorld_log =  fWorldPhysVol->GetLogicalVolume();
    G4int num_volumes = 9;
    auto *inner_volumes = new G4String[num_volumes];
    inner_volumes[0] = "element_log"; inner_volumes[1] = "Detector_log"; inner_volumes[2] = "Bipxl_log"; inner_volumes[3] = "membranepxl_log";
    inner_volumes[4] = "gridpiece_log"; inner_volumes[5] = "wafer_log"; inner_volumes[6] = "ACDpxl_log"; inner_volumes[7] = "AlFilter_log";
    inner_volumes[8] = "ACIC_log";

    G4cout << "Adding volumes to the InnerRegion" << G4endl;
    // Add volumes to the InnerRegion
    for (G4int j = 0; j <= num_volumes; j++)
    {
      for (G4int i = 0; i < (G4int)fWorld_log->GetNoDaughters(); i++)
      {
        if (fWorld_log->GetDaughter(i)->GetLogicalVolume()->GetName() == inner_volumes[j])
        {
          InnerRegion->AddRootLogicalVolume(fWorld_log->GetDaughter(i)->GetLogicalVolume());
          G4cout << "Added: " << fWorld_log->GetDaughter(i)->GetLogicalVolume()->GetName() << G4endl;
        }
      }
    }

    // Read volume names
    auto *pvs = G4PhysicalVolumeStore::GetInstance();
    if (pvs == nullptr)
    {
      G4cout << "PhysicalVolumeStore not accessible" << G4endl;
    }
    else
    {
      G4cout << "Volumes imported in the PhysicalVolumeStore. " << G4endl;
      G4int length = 0;
      length = pvs->size();
      for (G4int i= 0; i<length; i++)
      {
        G4String volName = ((*pvs)[i])->GetName();
      }
    }
    fWorld_phys = fWorldPhysVol;
  }

  // Set production cuts
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(200*eV, 1000*GeV);
  auto* InnerCuts = new G4ProductionCuts;
  InnerCuts->SetProductionCut(1*mm);
  InnerRegion->SetProductionCuts(InnerCuts);

  return fWorld_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimitAlfilta != nullptr)&&(maxStep>0.)){ fStepLimitAlfilta->SetMaxAllowedStep(maxStep); }
  if ((fStepLimitPolyfilta  != nullptr)&&(maxStep>0.)){ fStepLimitPolyfilta->SetMaxAllowedStep(maxStep); }
  if ((fStepLimitmembrane  != nullptr)&&(maxStep>0.)){ fStepLimitmembrane->SetMaxAllowedStep(maxStep); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetDetectorConstruction::SetReadFile(G4String input_geometry)
{
  fReadFile = std::move(input_geometry);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
