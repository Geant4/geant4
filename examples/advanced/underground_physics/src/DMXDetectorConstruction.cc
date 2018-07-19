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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"

#include "DMXScintSD.hh"
#include "DMXPmtSD.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()  
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DMXDetectorMessenger(this);

  theUserLimitsForRoom     = 0; 
  theUserLimitsForDetector = 0; 
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume
  
  //Zero the G4Cache objects to contain logical volumes
  LXeSD.Put(0);
  pmtSD.Put(0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction() 
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials() 
{

#include "DMXDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;
  

  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;



  // Universe - room wall - CONCRETE ****************************************

  //NB: measured INSIDE of lab, therefore have to add twice wall thickness
  G4double wallThick   = 24.*cm;
  G4double worldWidth  = 470.0*cm + 2.*wallThick; // "x"
  G4double worldLength = 690.0*cm + 2.*wallThick; // "y"
  G4double worldHeight = 280.0*cm + 2.*wallThick; // "z"

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, world_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "world_phys", world_log, NULL, false,0);

  //  G4VisAttributes* world_vat= new G4VisAttributes(white);
  world_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //world_vat->SetVisibility(true);
  //world_vat->SetVisibility(false);
  //world_log->SetVisAttributes(world_vat);


  // Lab Space - AIR ********************************************************

  G4double labWidth  = worldWidth  - 2.*wallThick; //X
  G4double labLength = worldLength - 2.*wallThick; //Y
  G4double labHeight = worldHeight - 2.*wallThick; //Z

  G4Box* lab_box = new G4Box
     ("lab_box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight );
  lab_log  = new G4LogicalVolume(lab_box, lab_mat, "lab_log");
  lab_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lab_phys", 
     lab_log, world_phys, false,0);

  G4VisAttributes* lab_vat= new G4VisAttributes(white);
  //  lab_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  lab_vat->SetVisibility(true);
  lab_vat->SetVisibility(false);
  lab_log->SetVisAttributes(lab_vat);

// include room furniture: **************************************************

#include "DMXDetectorRoom.icc"

  // Now start with detector assembly:

  // first LN2 cooling container: *******************************************

  G4double PosZ = -25.3*cm; // extra z-pos to correspond with height in lab

  G4double LN2jacketRadius    = 107.5*mm;
  G4double LN2jacketHeight    = 590.0*mm;
  G4double jacketHeight       = 680.0*mm;
  G4double jacketflangeHeight = 53.0*mm;
  G4double LN2PosZ            = 0.5*jacketHeight + 0.5*LN2jacketHeight 
                                + jacketflangeHeight + PosZ;

  G4Tubs* LN2jacket_tube = new G4Tubs("LN2jacket_tube",
     0.*cm, LN2jacketRadius, 0.5*LN2jacketHeight, 0.*deg, 360.*deg);
  LN2jacket_log  = new G4LogicalVolume
    (LN2jacket_tube, LN2jacket_mat, "LN2jacket_log");
  LN2jacket_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,LN2PosZ),
     "LN2jacket_phys", LN2jacket_log, lab_phys, false,0);

  G4VisAttributes* LN2jacket_vat = new G4VisAttributes(lgrey);
  // LN2jacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  // LN2jacket_vat->SetVisibility(true);
  LN2jacket_log->SetVisAttributes(LN2jacket_vat);

  // LN2jacket vacuum: **********************

  G4double LN2jacketMetalThick = 2.0*mm;
  G4double LN2vacuumRadius     = LN2jacketRadius - LN2jacketMetalThick;
  G4double LN2vacuumHeight     = LN2jacketHeight - LN2jacketMetalThick;

  G4Tubs* LN2vacuum_tube = new G4Tubs("LN2vacuum_tube",
     0.*cm, LN2vacuumRadius, 0.5*LN2vacuumHeight, 0.*deg, 360.*deg);
  LN2vacuum_log  = new G4LogicalVolume
    (LN2vacuum_tube, vacuum_mat, "LN2vacuum_log");
  LN2vacuum_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "LN2vacuum_phys", LN2vacuum_log, LN2jacket_phys, false,0);

  LN2vacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  // LN2 vessel: ************************************************************

  G4double LN2Radius       = 76.0*mm;
  G4double LN2Height       = 590.0*mm - 2.*LN2jacketMetalThick;
  G4double LN2vesselRadius = LN2Radius + LN2jacketMetalThick;
  G4double LN2vesselHeight = LN2Height;
  G4double LN2vesselPosZ   = 0.5*LN2vacuumHeight - 0.5*LN2vesselHeight;

  G4Tubs* LN2vessel_tube = new G4Tubs("LN2vessel_tube",
     0.*cm, LN2vesselRadius, 0.5*LN2vesselHeight, 0.*deg, 360.*deg);
  LN2vessel_log  = new G4LogicalVolume
    (LN2vessel_tube, LN2jacket_mat, "LN2vessel_log");
  LN2vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,LN2vesselPosZ),
     "LN2vessel_phys", LN2vessel_log, LN2vacuum_phys, false,0);

  G4VisAttributes* LN2vessel_vat = new G4VisAttributes(lgrey);
  // LN2vessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  // LN2vessel_vat->SetVisibility(true);
  LN2vessel_log->SetVisAttributes(LN2vessel_vat);


  // and finally LN2: *******************************************************

  G4Tubs* LN2_tube = new G4Tubs("LN2_tube",
     0.*cm, LN2Radius, 0.5*LN2Height, 0.*deg, 360.*deg);
  LN2_log  = new G4LogicalVolume(LN2_tube, LN2_mat, "LN2_log");
  LN2_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "LN2_phys", LN2_log, LN2vessel_phys, false,0);

  G4VisAttributes* LN2_vat = new G4VisAttributes(green);
  LN2_vat->SetVisibility(true);
  LN2_log->SetVisAttributes(LN2_vat);


  // outer vacuum jacket volume: stainless steel ****************************

  G4double jacketRadius     = 127.5*mm;
  //  G4double jacketHeight     = 680.0*mm; // defined above to get full-height
  G4double jacketMetalThick = LN2jacketMetalThick;

  G4Tubs* jacket_tube = new G4Tubs("jacket_tube",
     0.*cm, jacketRadius, 0.5*jacketHeight, 0.*deg, 360.*deg);
  jacket_log  = new G4LogicalVolume(jacket_tube, jacket_mat, "jacket_log");
  jacket_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,PosZ),
     "jacket_phys", jacket_log, lab_phys, false,0);

  G4VisAttributes* jacket_vat = new G4VisAttributes(grey);
  // jacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  jacket_log->SetVisAttributes(jacket_vat);


  // outer vacuum jacket flanges: stainless steel *************************

  G4double jacketflangeRadius  = 127.5*mm;
  // G4double jacketflangeHeight = 53.0*mm; // defined above to get full-height
  G4double topjacketflangePosZ = 0.5*(jacketHeight+jacketflangeHeight);

  G4Tubs* jacketflange_tube = new G4Tubs("jacketflange_tube",
     0.*cm, jacketflangeRadius, 0.5*jacketflangeHeight, 0.*deg, 360.*deg);
  jacketflange_log     = new G4LogicalVolume
    (jacketflange_tube, jacketflange_mat, "jacketflange_log");
  topjacketflange_phys = new G4PVPlacement
    (0, G4ThreeVector(0.,0.,topjacketflangePosZ + PosZ),
     "topjacketflange_phys", jacketflange_log, lab_phys, false,0);
  bottomjacketflange_phys = new G4PVPlacement
    (0, G4ThreeVector(0.,0.,-topjacketflangePosZ + PosZ),
     "bottomjacketflange_phys", jacketflange_log, lab_phys, false,0);

  // jacketflange_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  jacketflange_log->SetVisAttributes(jacket_vat);

  // vacuum **************************************************************

  G4double vacuumRadius = jacketRadius - jacketMetalThick;
  G4double vacuumHeight = jacketHeight - jacketMetalThick;

  G4Tubs* vacuum_tube = new G4Tubs("vacuum_tube",
     0.*cm, vacuumRadius, 0.5*vacuumHeight, 0.*deg, 360.*deg);
  vacuum_log  = new G4LogicalVolume(vacuum_tube, vacuum_mat, "vacuum_log");
  vacuum_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "vacuum_phys", vacuum_log, jacket_phys, false,0);

  // G4VisAttributes* vacuum_vat= new G4VisAttributes(lgrey);
  vacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());


  // copper cooling jacket volume: **************************************

  G4double copperMetalThick = 3.0*mm;
  G4double copperRadius     = 103.5*mm + copperMetalThick;
  G4double copperHeight     = 420.0*mm;
  G4double copperInner      = copperRadius - copperMetalThick;
  G4double vesselHeight     = 320.0*mm;
  G4double copperVPos       = 0.5*(vesselHeight-copperHeight) + 13.0*cm;
  G4double coppertopThick   = 1.0*cm;
  G4double coppertopVPos    = copperVPos + 0.5*(coppertopThick+copperHeight);

  G4Tubs* copper_tube = new G4Tubs("copper_tube",
     copperInner, copperRadius, 0.5*copperHeight, 0.*deg, 360.*deg);
  copper_log  = new G4LogicalVolume(copper_tube, copper_mat, "copper_log");
  copper_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,copperVPos), 
     "copper_phys", copper_log, vacuum_phys, false,0);

  G4Tubs* coppertop_tube = new G4Tubs("coppertop_tube",
     0.*cm, copperRadius, 0.5*coppertopThick, 0.*deg, 360.*deg);
  coppertop_log  = new G4LogicalVolume
    (coppertop_tube, copper_mat, "coppertop_log");  
  coppertop_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,coppertopVPos), 
     "coppertop_phys", coppertop_log, vacuum_phys, false,0);

  G4VisAttributes* copper_vat = new G4VisAttributes(orange);
  //  copper_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  copper_log->SetVisAttributes(copper_vat);
  coppertop_log->SetVisAttributes(copper_vat);

  // inner vessel jacket volume: stainless steel ************************

  //  G4double vesselHeight = 320.0*mm; // - moved earlier
  G4double vesselMetalThick      = jacketMetalThick;
  G4double vesselRadius          = 75.0*mm + vesselMetalThick;
  G4double vesselflangeRadius    = 101.5*mm;
  G4double vesselflangeThick     = 40.0*mm;
  G4double PMTvesselRadius       = 31.0*mm + vesselMetalThick;
  G4double PMTvesselHeight       = 152.0*mm;
  G4double pmtvesselflangeRadius = 52.0*mm;
  G4double pmtvesselflangeThick  = 32.0*mm;
  G4double vesselVPos            = 7.0*cm;
  G4double TotalvesselHeight     = PMTvesselHeight + vesselHeight;

  G4Tubs* vessel_tube    = new G4Tubs("vessel_tube",
     0.*cm, vesselRadius, 0.5*vesselHeight, 0.*deg, 360.*deg);
  G4Tubs* PMTvessel_tube = new G4Tubs("PMTvessel_tube",
     0.*cm, PMTvesselRadius, 0.5*PMTvesselHeight, 0.*deg, 360.*deg);

  G4UnionSolid* vessel_sol = new G4UnionSolid
    ("vessel_sol", vessel_tube, PMTvessel_tube,
     G4Transform3D(G4RotationMatrix(), 
		   G4ThreeVector(0,0,-0.5*(vesselHeight+PMTvesselHeight))));

  vessel_log  = new G4LogicalVolume(vessel_sol, vessel_mat, "vessel_log");
  vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,vesselVPos), 
     "vessel_phys", vessel_log, vacuum_phys, false,0);


  // flanges: 1=upper half (diff. inner diam.) 2=lower half
  G4Tubs* vesseltop_flange1 = new G4Tubs("vesseltop_flange1",
     0.*cm, vesselflangeRadius, 0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesseltop_log1  = new G4LogicalVolume
    (vesseltop_flange1, vessel_mat, "vesseltop_log1");  
  vesseltop_phys1 = new G4PVPlacement
    (0, 
     G4ThreeVector(0.,0.,0.5*(vesselHeight+0.5*vesselflangeThick)+vesselVPos), 
     "vesseltop_phys1", vesseltop_log1, vacuum_phys, false,0);

  G4Tubs* vesseltop_flange2 = new G4Tubs("vesseltop_flange2",vesselRadius, 
    vesselflangeRadius, 0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesseltop_log2  = new G4LogicalVolume
    (vesseltop_flange2, vessel_mat, "vesseltop_log2");  
  vesseltop_phys2 = new G4PVPlacement
    (0, 
     G4ThreeVector(0.,0.,0.5*(vesselHeight-0.5*vesselflangeThick)+vesselVPos), 
     "vesseltop_phys2", vesseltop_log2, vacuum_phys, false,0);


  G4Tubs* vesselbottom_flange1 = new G4Tubs
    ("vesselbottom_flange1",vesselRadius, vesselflangeRadius, 
     0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesselbottom_log1  = new G4LogicalVolume
    (vesselbottom_flange1, vessel_mat, "vesselbottom_log1");  
  vesselbottom_phys1 = new G4PVPlacement(0, 
     G4ThreeVector(0.,0.,-0.5*(vesselHeight-0.5*vesselflangeThick)+vesselVPos),
     "vesselbottom_phys1", vesselbottom_log1, vacuum_phys, false,0);

  G4Tubs* vesselbottom_flange2 = new G4Tubs
    ("vesselbottom_flange2",PMTvesselRadius, vesselflangeRadius, 
     0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesselbottom_log2  = new G4LogicalVolume
    (vesselbottom_flange2, vessel_mat, "vesselbottom_log2");  
  vesselbottom_phys2 = new G4PVPlacement(0, 
     G4ThreeVector(0.,0.,-0.5*(vesselHeight+0.5*vesselflangeThick)+vesselVPos),
     "vesselbottom_phys2", vesselbottom_log2, vacuum_phys, false,0);


  G4Tubs* pmtvesselbottom_flange1 = new G4Tubs
    ("pmtvesselbottom_flange1", PMTvesselRadius, pmtvesselflangeRadius, 
     0.25*pmtvesselflangeThick, 0.*deg, 360.*deg);
  pmtvesselbottom_log1  = new G4LogicalVolume
    (pmtvesselbottom_flange1, vessel_mat, "pmtvesselbottom_log1");  
  pmtvesselbottom_phys1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,
    (-0.5*vesselHeight-PMTvesselHeight+vesselVPos+0.25*pmtvesselflangeThick)),
     "pmtvesselbottom_phys1", pmtvesselbottom_log1, vacuum_phys, false,0);

  G4Tubs* pmtvesselbottom_flange2 = new G4Tubs
    ("pmtvesselbottom_flange2", 0.*cm, pmtvesselflangeRadius, 
     0.25*pmtvesselflangeThick, 0.*deg, 360.*deg);
  pmtvesselbottom_log2  = new G4LogicalVolume
    (pmtvesselbottom_flange2, vessel_mat, "pmtvesselbottom_log2");  
  pmtvesselbottom_phys2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,
     -0.5*vesselHeight-PMTvesselHeight+vesselVPos-0.25*pmtvesselflangeThick),
     "pmtvesselbottom_phys2", pmtvesselbottom_log2, vacuum_phys, false,0);


  G4VisAttributes* vessel_vat     = new G4VisAttributes(grey);
  G4VisAttributes* pmtvessel_vat  = new G4VisAttributes(yellow);
  G4VisAttributes* pmtvessel_vat2 = new G4VisAttributes(green);
  //  vessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  vessel_vat->SetForceSolid(true);
  //  pmtvessel_vat->SetForceSolid(true);
  //  pmtvessel_vat2->SetForceSolid(true);
  vessel_log->SetVisAttributes(vessel_vat);
  vesseltop_log1->SetVisAttributes(vessel_vat);
  vesselbottom_log1->SetVisAttributes(vessel_vat);
  vesseltop_log2->SetVisAttributes(pmtvessel_vat);
  vesselbottom_log2->SetVisAttributes(pmtvessel_vat);
  //  pmtvesselbottom_log->SetVisAttributes(vessel_vat);
  pmtvesselbottom_log1->SetVisAttributes(vessel_vat);
  pmtvesselbottom_log2->SetVisAttributes(pmtvessel_vat2);



  // *********************************************************************
  // grid#1 to mirror surface: 21.75 mm
  // LXe height = 15.75 mm, gXe height = 6.00 mm
  // NB: Increased liquid height by 1mm - to take away problem with 
  // over-lapping volumes/ring pronounced from liquid phase..........
  // *********************************************************************

  // detector volume: gas phase ******************************************

  G4double mirrorVPos     = 21.3*cm;
  G4double gasGap         = 6.0*mm;
  G4double DetectorRadius = vesselRadius - vesselMetalThick;
  G4double GXeHeight      = TotalvesselHeight - mirrorVPos + gasGap;
  G4double GXeVPos        = 0.5*vesselHeight - 0.5*GXeHeight;

  G4Tubs* GXe_tube = new G4Tubs("GXe_tube",
     0.*cm, DetectorRadius, 0.5*GXeHeight, 0.*deg, 360.*deg);
  GXe_log  = new G4LogicalVolume(GXe_tube, GXe_mat, "GXe_log");
  GXe_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,GXeVPos), 
     "GXe_phys", GXe_log, vessel_phys, false,0);

  G4VisAttributes* GXe_vat = new G4VisAttributes(cyan);
  // GXe_vat->SetForceSolid(true);
  GXe_vat->SetVisibility(true);
  GXe_log->SetVisAttributes(GXe_vat);


  // liquid phase *******************************************************

  G4double LXeHeight         = mirrorVPos - gasGap;
  G4double PMTDetectorRadius = PMTvesselRadius - vesselMetalThick;
  G4double PMTDetectorHeight = PMTvesselHeight;
  G4double LXeTubeHeight     = LXeHeight - PMTDetectorHeight;
  G4double LXe_solVPos       = -0.5*(LXeTubeHeight+PMTDetectorHeight);
  G4double LXeVPos           = -0.5*TotalvesselHeight + 0.5*LXeHeight;

  G4Tubs* LXe_tube = new G4Tubs("GXe_tube",
     0.*cm, DetectorRadius, 0.5*LXeTubeHeight, 0.*deg, 360.*deg);
  G4Tubs* PMTdetector_tube = new G4Tubs("PMTdetector_tube",
   0.*cm, PMTDetectorRadius, 0.5*PMTDetectorHeight, 0.*deg, 360.*deg);

  G4UnionSolid* LXe_sol = new G4UnionSolid
    ("LXe_sol", LXe_tube, PMTdetector_tube,
    G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,LXe_solVPos)));

  LXe_log  = new G4LogicalVolume(LXe_sol, LXe_mat, "LXe_log");
  LXe_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, LXeVPos), 
    "LXe_phys", LXe_log, vessel_phys, false, 0);

  // attributes
  G4VisAttributes* LXe_vat = new G4VisAttributes(lblue);
  // LXe_vat->SetForceSolid(true);
  LXe_vat->SetVisibility(true);
  LXe_log->SetVisAttributes(LXe_vat);

  
  // Gas phase vessel lagging - for optical properties:

  G4double laggingThickness = 10.*micrometer;
  G4double laggingRadius    = DetectorRadius - laggingThickness;

  G4Tubs* gaslag_tube = new G4Tubs("gaslag_tube", laggingRadius, 
     DetectorRadius, 0.5*GXeHeight, 0.*deg, 360.*deg);
  gaslag_log  = new G4LogicalVolume(gaslag_tube, vessel_mat, "gaslag_log");
  gaslag_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm),
    "gaslag_phys", gaslag_log, GXe_phys, false, 0);

  // attributes
  G4VisAttributes* gaslag_vat = new G4VisAttributes(lgreen);
  // gaslag_vat->SetForceSolid(true);
  gaslag_vat->SetVisibility(true);
  gaslag_log->SetVisAttributes(gaslag_vat);


  // liquid phase vessel lagging - for optical properties:

  G4double lagTubeRadius = DetectorRadius - laggingThickness;
  G4double lagTubeHeight = LXeHeight - PMTDetectorHeight;
  G4double lagPMTRadius  = PMTDetectorRadius - laggingThickness;
  G4double lagPMTHeight  = PMTDetectorHeight;

  G4Tubs* liqLag_tube = new G4Tubs("liqlag_tube", lagTubeRadius,
     DetectorRadius, 0.5*lagTubeHeight, 0.*deg, 360.*deg);
  G4Tubs* lagPMT_tube = new G4Tubs("lagPMT_tube", lagPMTRadius, 
     PMTDetectorRadius, 0.5*lagPMTHeight, 0.*deg, 360.*deg);

  G4UnionSolid* liqLag_sol = new G4UnionSolid
    ("liqLag_sol", liqLag_tube, lagPMT_tube,
    G4Transform3D(G4RotationMatrix(),G4ThreeVector(0,0,LXe_solVPos)));

  liqLag_log  = new G4LogicalVolume(liqLag_sol, vessel_mat, "liqLag_log");
  liqLag_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm), 
    "liqLag_phys", liqLag_log, LXe_phys, false, 0);

  // attributes
  G4VisAttributes* liqLag_vat = new G4VisAttributes(magenta);
  // liqLag_vat->SetForceSolid(true);
  liqLag_vat->SetVisibility(true);
  liqLag_log->SetVisAttributes(liqLag_vat);


  // Vessel Wall Optical Surface definition:
  G4OpticalSurface* OpVesselSurface = new G4OpticalSurface
    ("VesselSurface", unified, polished, dielectric_metal);

  // created optical lagging onto vessel - to avoid clash between over-lapping
  // liquid and gas phase - so removed below:
  /*
  G4LogicalBorderSurface* VesselSurface;
  VesselSurface = new G4LogicalBorderSurface
    ("Vessel", liqPhase_phys, vessel_phys, OpVesselSurface);
  */

  const G4int NUM = 2;
  G4double vessel_PP[NUM]   = { 6.5*eV, 7.50*eV };
  G4double vessel_REFL[NUM] = { 0.2, 0.2 };
  G4MaterialPropertiesTable* vessel_mt = new G4MaterialPropertiesTable();
  vessel_mt->AddProperty("REFLECTIVITY", vessel_PP, vessel_REFL, NUM);
  OpVesselSurface->SetMaterialPropertiesTable(vessel_mt);

  // G4LogicalBorderSurface* VesselTopSurface = 
  new G4LogicalBorderSurface
    ("VesselTop", GXe_phys, vesseltop_phys1, OpVesselSurface);

  //G4LogicalBorderSurface* VesselBottomSurface = 
  new G4LogicalBorderSurface
    ("VesselBottom", LXe_phys, vesselbottom_phys2, OpVesselSurface);

  //G4LogicalBorderSurface* GasVesselSurface = 
  new G4LogicalBorderSurface
    ("GasVessel", GXe_phys, gaslag_phys, OpVesselSurface);

  //G4LogicalBorderSurface* LiquidVesselSurface =
  new G4LogicalBorderSurface
    ("LiquidVessel", LXe_phys, liqLag_phys, OpVesselSurface);



  // Cu Shield **********************************************************

  G4double CuShieldHeight      = 17.7*cm;
  G4double CuShieldThickness   = 2.4*mm;
  G4double CuShieldOuterRadius = 3.0*cm;
  G4double CuShieldInnerRadius = CuShieldOuterRadius-CuShieldThickness;
  G4double CuShieldVPosition   = -0.5*LXeTubeHeight - PMTDetectorHeight 
                                + 0.5*CuShieldHeight;

  // Zero co-ordinate of the union is the zero of the first volume, 
  // i.e. the offset is still present

  G4Tubs* CuShield_tube = new G4Tubs("CuShield_tube", CuShieldInnerRadius,
     CuShieldOuterRadius, 0.5*CuShieldHeight, 0.*deg, 360.*deg);
  CuShield_log  = new G4LogicalVolume(CuShield_tube, CuShield_mat, 
				     "CuShield_log");
  CuShield_phys = new G4PVPlacement(0, 
     G4ThreeVector(0.*cm, 0.*cm, CuShieldVPosition), 
     "CuShield_phys", CuShield_log, LXe_phys, false, 0);

  //  G4VisAttributes* CuShield_vat= new G4VisAttributes(magenta);
  G4VisAttributes* CuShield_vat = new G4VisAttributes(brown);
  //  CuShield_vat->SetForceSolid(true);
  CuShield_vat->SetVisibility(true);
  CuShield_log->SetVisAttributes(CuShield_vat);

  // Cu shield surface
  G4double sigalpha;
  G4OpticalSurface* OpCuShieldSurface = new G4OpticalSurface
    ("ShieldSurface", unified, ground, dielectric_metal, sigalpha=30.0*deg);
  //G4LogicalBorderSurface* ShieldSurface = 
  new G4LogicalBorderSurface
    ("Shield", LXe_phys, CuShield_phys, OpCuShieldSurface);

  G4double CuShield_PP[NUM]   = { 7.0*eV, 7.50*eV };
  G4double CuShield_REFL[NUM] = { 0.3, 0.2 };
  G4MaterialPropertiesTable *CuShield_mt = new G4MaterialPropertiesTable();
  CuShield_mt->AddProperty("REFLECTIVITY", CuShield_PP, CuShield_REFL, NUM);
  OpCuShieldSurface->SetMaterialPropertiesTable(CuShield_mt);


  // rings ***************************************************************

  G4double ringHeight      =  4.*mm;
  G4double ringOuterRadius =  4.0*cm;
  G4double ringInnerRadius =  CuShieldOuterRadius;
  G4double ringVOffset     =  0.5*ringHeight;
  G4double ringVPosition   =  -0.5*GXeHeight + gasGap +ringVOffset;

  G4Tubs* ring_tube=new G4Tubs("ring_tube", ringInnerRadius,
     ringOuterRadius, 0.5*ringHeight, 0.*deg, 360.*deg);
  ring_log = new G4LogicalVolume(ring_tube, ring_mat, "ring_log");

  // optical surface: ring materials table
  G4double ring_PP[NUM]   = { 6.00*eV, 7.50*eV };
  G4double ring_REFL[NUM] = { 0.7, 0.65 };
  G4MaterialPropertiesTable *ring_mt = new G4MaterialPropertiesTable();
  ring_mt->AddProperty("REFLECTIVITY", ring_PP, ring_REFL, NUM);

  G4OpticalSurface* OpRingSurface = new G4OpticalSurface
    ("RingSurface", unified, ground, dielectric_metal, sigalpha=10.*deg);
  // last argument is surface roughness if it's non-polished - i.e. ground
  OpRingSurface->SetMaterialPropertiesTable(ring_mt);

  // rings inside gas phase
  ring_phys_gas[0] = new G4PVPlacement(0, G4ThreeVector
    (0.*cm, 0.*cm, ringVPosition),"ring_phys0",ring_log,GXe_phys,false, 0);
  //G4LogicalBorderSurface* RingSurface_gas0 = 
  new G4LogicalBorderSurface
    ("Ring", GXe_phys, ring_phys_gas[0], OpRingSurface);

  ring_phys_gas[1] = new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.0*mm),
     "ring_phys1",ring_log, GXe_phys, false, 0);
  //G4LogicalBorderSurface* RingSurface_gas1 = 
  new G4LogicalBorderSurface
     ("Ring", GXe_phys, ring_phys_gas[1], OpRingSurface);


  // rings inside liquid phase:
  ringVPosition = 0.5*LXeTubeHeight;

  ring_phys_liq[0] = new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=0.5*ringHeight),
     "ring_phys2",ring_log,LXe_phys, false, 0);
  //G4LogicalBorderSurface* RingSurface_liq0 = 
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[0], OpRingSurface);

  ring_phys_liq[1] = new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys3",ring_log, LXe_phys, false, 0);
  //G4LogicalBorderSurface* RingSurface_liq1 =
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[1], OpRingSurface);

  ring_phys_liq[2]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight),
     "ring_phys4",ring_log, LXe_phys, false, 0);
  //G4LogicalBorderSurface* RingSurface_liq2 =
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[2], OpRingSurface);

  ring_phys_liq[3]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight),
     "ring_phys5",ring_log, LXe_phys, false, 0);
  //G4LogicalBorderSurface* RingSurface_liq3 =
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[3], OpRingSurface);

  ring_phys_liq[4]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys6",ring_log, LXe_phys,false, 0);
  //G4LogicalBorderSurface* RingSurface_liq4 =
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[4], OpRingSurface);

  ring_phys_liq[5]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys7",ring_log, LXe_phys,false, 0);
  //G4LogicalBorderSurface* RingSurface_liq5 =
  new G4LogicalBorderSurface
    ("Ring", LXe_phys, ring_phys_liq[5], OpRingSurface);


  G4VisAttributes* ring_vat= new G4VisAttributes(lgrey);
  ring_vat->SetVisibility(true);
  ring_log->SetVisAttributes(ring_vat);


  // Mirror *************************************************************

  G4double mirrorHeight    = 2.0*mm;
  G4double mirrorRadius    = ringInnerRadius;
  G4double mirrorVOffset   = 0.5*ringHeight;
  G4double mirrorVPosition = -0.5*GXeHeight + gasGap +mirrorVOffset;

  G4Tubs* mirror_tube = new G4Tubs("mirror_tube", 0.*cm, mirrorRadius,
     0.5*mirrorHeight, 0.*deg, 360.*deg);
  mirror_log  = new G4LogicalVolume(mirror_tube, mirror_mat, "mirror_log");
  mirror_phys = new G4PVPlacement(0, 
     G4ThreeVector(0.*cm, 0.*cm, mirrorVPosition),
     "mirror_phys", mirror_log, GXe_phys, false, 0);

  G4VisAttributes* mirror_vat = new G4VisAttributes(cyan);
  mirror_vat->SetVisibility(true);
  //  mirror_vat->SetForceSolid(true);
  mirror_log->SetVisAttributes(mirror_vat);


  // mirror surface
  G4OpticalSurface * OpMirrorSurface = new G4OpticalSurface
    ("MirrorSurface", unified, ground, dielectric_metal, sigalpha=5.0*deg);
  //G4LogicalBorderSurface* MirrorSurface = 
  new G4LogicalBorderSurface
    ("Mirror", GXe_phys, mirror_phys, OpMirrorSurface);

  G4double mirror_PP[NUM]   = { 6.00*eV, 7.50*eV };
  G4double mirror_REFL[NUM] = { 0.83, 0.78 };
  G4MaterialPropertiesTable *mirror_mt = new G4MaterialPropertiesTable();
  mirror_mt->AddProperty("REFLECTIVITY", mirror_PP, mirror_REFL, NUM);
  OpMirrorSurface->SetMaterialPropertiesTable(mirror_mt);


  // Grids  *************************************************************

  G4double gridHeight     = 0.100*mm;
  G4double gridRadius     = ringInnerRadius;
  G4double grid1VOffset   = 3.5*ringHeight+1.75*mm;
  G4double grid1VPosition = 0.5*LXeTubeHeight - grid1VOffset;
  G4double grid2VOffset   = 4.5*ringHeight+3.50*mm;
  G4double grid2VPosition = 0.5*LXeTubeHeight - grid2VOffset;

  G4Tubs* grid_tube = new G4Tubs("grid_tube", 0.*cm, gridRadius,
     0.5*gridHeight, 0.*deg, 360.*deg);

  grid1_log  = new G4LogicalVolume(grid_tube, grid_mat, "grid1_log");
  grid1_phys = new G4PVPlacement(0,G4ThreeVector(0.*cm, 0.*cm, grid1VPosition),
     "grid1_phys", grid1_log, LXe_phys, false, 0);
  grid2_log  = new G4LogicalVolume(grid_tube, grid_mat, "grid2_log");
  grid2_phys = new G4PVPlacement(0,G4ThreeVector(0.*cm, 0.*cm, grid2VPosition),
     "grid2_phys", grid2_log, LXe_phys, false, 0);

  G4VisAttributes* grid_vat = new G4VisAttributes(red);
  grid_vat->SetVisibility(true);
  grid1_log->SetVisAttributes(grid_vat);
  grid2_log->SetVisAttributes(grid_vat);
  
  
  // alpha source holder ************************************************
  
  G4double alphaHeight     = 0.7*mm; // total lead thickness = 1.23 mm
  G4double recessHeight    = 0.3*mm;  // totals lead thickness = 1.23 mm
  G4double alphaRadius     = 1.2*mm; // was 0.8
  G4double recessRadius    = 0.35*mm; // was 0.5 was 0.2
  G4double recessVPosition = 0.5*(alphaHeight - recessHeight);

  G4double alphaVOffset    = grid1VOffset-0.5*alphaHeight - gridHeight;
  G4double americiumHeight = 3000.*nanometer; // assumes ~1% Am  
  G4double extra_offset    = (recessHeight +0.3*mm + 0.5*americiumHeight); 
  G4double alphaVPosition  = 0.5*LXeTubeHeight - alphaVOffset + extra_offset;

  G4Tubs* alpha_tube  = new G4Tubs("alpha_tube", 0.*cm, alphaRadius,
     0.5*alphaHeight,  0.*deg, 360.*deg);
  G4Tubs* recess_tube = new G4Tubs("recess_tube", 0.*cm, recessRadius,
     0.5*recessHeight, 0.*deg, 360.*deg);

  G4SubtractionSolid* alpha_sol = new G4SubtractionSolid
    ("alpha_sol", alpha_tube, recess_tube, G4Transform3D
     (G4RotationMatrix(), G4ThreeVector(0. ,0. , recessVPosition)));
  alpha_log  = new G4LogicalVolume(alpha_sol, alpha_mat, "alpha_log");

  alpha_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., alphaVPosition),
                         "alpha_phys", alpha_log, LXe_phys, false, 0);

  G4VisAttributes* alpha_vat = new G4VisAttributes(white);
  alpha_vat->SetVisibility(true);
  alpha_log ->SetVisAttributes(alpha_vat);

  // alpha source HOLDER surface
  G4OpticalSurface* OpAlphaSurface = new G4OpticalSurface("AlphaSurface", 
  unified, ground, dielectric_metal, sigalpha=20.0*deg);
  //G4LogicalBorderSurface* AlphaSurface =
  new G4LogicalBorderSurface
    ("Alpha", LXe_phys, alpha_phys, OpAlphaSurface);

  G4double alpha_PP[NUM]   = { 6.00*eV, 7.50*eV };
  G4double alpha_REFL[NUM] = { 0.05, 0.05 };
  G4MaterialPropertiesTable *alpha_mt = new G4MaterialPropertiesTable();
  alpha_mt->AddProperty("REFLECTIVITY", alpha_PP, alpha_REFL, NUM);
  OpAlphaSurface->SetMaterialPropertiesTable(alpha_mt);

  // americium ***********************************************************

  // moved above for the "extra_offset":
  // G4double americiumHeight    = 600.*nanometer; // assumes ~1% Am
  G4double americiumRadius    = recessRadius - 50.0*micrometer;
  G4double americiumVOffset   = 0.5*(alphaHeight-americiumHeight)-recessHeight;
  G4double americiumVPosition = americiumVOffset;

  sourceZ = vesselVPos + LXeVPos + alphaVPosition + americiumVPosition + PosZ;
  G4cout << G4endl << "Calibration source centre (X,Y,Z):  0, 0, " 
	 << sourceZ/mm << " mm" << G4endl;

  G4Tubs* americium_tube = new G4Tubs("americium_tube", 0.*cm,
     americiumRadius, 0.5*americiumHeight, 0.*deg, 360.*deg);
  americium_log  = new G4LogicalVolume(americium_tube, americium_mat,
     "americium_log");
  americium_phys = new G4PVPlacement(0, G4ThreeVector(0., 0.,
     americiumVPosition),"americium_phys", americium_log, alpha_phys,false,0);

  // americium optical properties:
  G4OpticalSurface* OpAmericiumSurface = new G4OpticalSurface
    ("AmericiumSurface", unified, ground, dielectric_metal, sigalpha=5.0*deg);
  //G4LogicalBorderSurface* AmericiumSurface =
  new G4LogicalBorderSurface
    ("Americium", LXe_phys, americium_phys, OpAmericiumSurface);

  G4double americium_PP[NUM]   = { 6.00*eV, 7.50*eV };
  G4double americium_REFL[NUM] = { 0.7, 0.65 };
  G4MaterialPropertiesTable *americium_mt = new G4MaterialPropertiesTable();
  americium_mt->AddProperty("REFLECTIVITY", americium_PP, americium_REFL, NUM);
  OpAlphaSurface->SetMaterialPropertiesTable(americium_mt);

  G4VisAttributes* americium_vat= new G4VisAttributes(cyan);
  americium_vat->SetVisibility(true);
  americium_vat->SetForceSolid(true);
  americium_log->SetVisAttributes(americium_vat);


  // Photomultiplier: ETL 9829 QA ****************************************

  G4double pmtHeight    = 12.0*cm;
  G4double pmtRadius    = 2.6*cm;
  G4double pmtVOffset   = 1.0*cm;
  G4double pmtVPosition = -0.5*(LXeTubeHeight+pmtHeight)+pmtVOffset;

  G4Sphere* pmt_window = new G4Sphere("pmt_sphere", 0.*cm, 2.*pmtRadius, 
     0.*deg, 360.*deg, 0.*deg, 30.0*deg);
  G4Tubs* pmt_tube = new G4Tubs("pmt_tube", 0.*cm,  pmtRadius, 0.5*pmtHeight,
     0.*deg, 360.*deg);
  
  G4UnionSolid* pmt_sol = new G4UnionSolid("pmt_sol", pmt_tube, pmt_window,
    G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0.5*pmtHeight
    -2.*pmtRadius*std::cos(30.0*deg))));

  pmt_log  = new G4LogicalVolume(pmt_sol, pmt_mat, "pmt_log");
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(0.*cm, 0.*cm, pmtVPosition),
     "pmt_phys", pmt_log, LXe_phys, false, 0);

  G4OpticalSurface* pmt_opsurf = new G4OpticalSurface
    ("pmt_opsurf",unified, polished, dielectric_dielectric);
  //G4LogicalBorderSurface* pmt_surf = 
  new G4LogicalBorderSurface
    ("pmt_surf", LXe_phys, pmt_phys, pmt_opsurf);

  G4VisAttributes* pmt_vat= new G4VisAttributes(blue);
  pmt_vat->SetForceSolid(true);
  pmt_vat->SetVisibility(true);
  pmt_log->SetVisAttributes(pmt_vat);


  // photocathode *******************************************************

  G4double phcathVOffset     = 0.5*pmtHeight-2.*pmtRadius*std::cos(30.0*deg);
  G4double phcathVPosition   = phcathVOffset;

  G4Sphere* phcath_sol = new G4Sphere("phcath_sphere",
     2.*pmtRadius-1.6*mm, 2.*pmtRadius-1.59*mm, 0.*deg, 360.*deg, 0.*deg, 
     27.0*deg);

  phcath_log  = new G4LogicalVolume(phcath_sol, phcath_mat, "phcath_log");
  phcath_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., phcathVPosition),
     "phcath_phys", phcath_log, pmt_phys, false, 0);

  G4OpticalSurface*  phcath_opsurf = new G4OpticalSurface("phcath_opsurf",
     unified, polished, dielectric_dielectric);
  //G4LogicalBorderSurface* phcath_surf = 
  new G4LogicalBorderSurface
    ("phcath_surf", pmt_phys, phcath_phys, phcath_opsurf);

  G4double phcath_PP[NUM]   = { 6.00*eV, 7.50*eV };
  // G4double phcath_REFL[NUM] = { 0.0, 0.0};
  // G4MaterialPropertiesTable* phcath_mt = new G4MaterialPropertiesTable();
  // phcath_mt->AddProperty("REFLECTIVITY", phcath_PP, phcath_REFL, NUM);
  // phcath_opsurf->SetMaterialPropertiesTable(phcath_mt);


  //**Photocathode surface properties
  G4double photocath_EFF[NUM]={1.,1.}; //Enables 'detection' of photons
  G4double photocath_ReR[NUM]={1.92,1.92};
  G4double photocath_ImR[NUM]={1.69,1.69};
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",phcath_PP,photocath_EFF,NUM);
  photocath_mt->AddProperty("REALRINDEX",phcath_PP,photocath_ReR,NUM);
  photocath_mt->AddProperty("IMAGINARYRINDEX",phcath_PP,photocath_ImR,NUM);
  G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);


  G4VisAttributes* phcath_vat= new G4VisAttributes(lblue);
  phcath_vat->SetForceSolid(true);
  phcath_vat->SetVisibility(true);
  phcath_log->SetVisAttributes(phcath_vat);

  new G4LogicalSkinSurface("photocath_surf",phcath_log,photocath_opsurf);

  // ......................................................................
  // attach user limits ...................................................

  
  G4cout << G4endl << "User Limits: " << G4endl 
	 << "\t theMaxTimeCuts:     " << G4BestUnit(theMaxTimeCuts,"Time")  
	 << G4endl
	 << "\t theRoomTimeCut:     " << G4BestUnit(theRoomTimeCut,"Time")  
	 << G4endl
	 << "\t theMaxStepSize:     " << G4BestUnit(theMaxStepSize,"Length")
	 << G4endl
	 << "\t theMinEKine:        " << G4BestUnit(theMinEkine,"Energy")   
	 << G4endl
	 << "\t minRoomMinEKine:    " << G4BestUnit(theRoomMinEkine,"Energy")
	 << G4endl << G4endl;

  if (theUserLimitsForRoom != 0) delete theUserLimitsForRoom;
  if (theUserLimitsForDetector != 0) delete theUserLimitsForDetector;

  theUserLimitsForRoom = new G4UserLimits(theMaxStepSize,   // step length max
					  DBL_MAX,          // track length max
					  theRoomTimeCut,   // Time cut
					  theRoomMinEkine); // min energy

#include "DMXDetectorRoomLimits.icc"

  theUserLimitsForDetector = new G4UserLimits(theDetectorStepSize,
					      DBL_MAX, // Track Max
					      theMaxTimeCuts,
					      theMinEkine);

      world_log->SetUserLimits(theUserLimitsForRoom);
        lab_log->SetUserLimits(theUserLimitsForRoom);
     jacket_log->SetUserLimits(theUserLimitsForRoom);
     vacuum_log->SetUserLimits(theUserLimitsForRoom);
     vessel_log->SetUserLimits(theUserLimitsForRoom);
        GXe_log->SetUserLimits(theUserLimitsForDetector);
	//        LXe_log->SetUserLimits(theUserLimitsForXenon);
        LXe_log->SetUserLimits(theUserLimitsForDetector);
   CuShield_log->SetUserLimits(theUserLimitsForDetector);
       ring_log->SetUserLimits(theUserLimitsForDetector);
     mirror_log->SetUserLimits(theUserLimitsForDetector);
      grid1_log->SetUserLimits(theUserLimitsForDetector);
      grid2_log->SetUserLimits(theUserLimitsForDetector);
      alpha_log->SetUserLimits(theUserLimitsForDetector);
  americium_log->SetUserLimits(theUserLimitsForDetector);
        pmt_log->SetUserLimits(theUserLimitsForDetector);
     phcath_log->SetUserLimits(theUserLimitsForDetector);

 return world_phys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
  // ......................................................................
  // sensitive detectors ..................................................
  // ......................................................................

  if (LXeSD.Get() == 0) 
    {    
      G4String name="/DMXDet/LXeSD";
      DMXScintSD* aSD = new DMXScintSD(name);
      LXeSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get());  
  if (LXe_log)    
    SetSensitiveDetector(LXe_log,LXeSD.Get());

  if (pmtSD.Get() == 0)
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get()); 
  if (phcath_log)
    SetSensitiveDetector(phcath_log,pmtSD.Get());

  return;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetRoomEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for ROOM
  theRoomMinEkine = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMinEkine(val); 
      G4cout << " Changing Room energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for Xenon Detector
  theMinEkine = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMinEkine(val);
      G4cout << "Changing Detector energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetRoomTimeCut(G4double val)
{
  // set room time cut:
  theRoomTimeCut = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMaxTime(val);
      G4cout << " Changing Room Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetTimeCut(G4double val)
{
  // set detector time cut:
  theMaxTimeCuts = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMaxTime(val);
      G4cout << " Changing Detector Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



