
#include "ExN04DetectorConstruction.hh"
#include "ExN04TrackerSD.hh"
#include "ExN04CalorimeterSD.hh"
#include "ExN04CalorimeterROGeometry.hh"
#include "ExN04MuonSD.hh"
#include "ExN04TrackerParametrisation.hh"
#include "ExN04CalorimeterParametrisation.hh"
#include "ExN04Field.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

ExN04DetectorConstruction::ExN04DetectorConstruction()
{

#include "ExN04DetectorParameterDef.icc"

}

ExN04DetectorConstruction::~ExN04DetectorConstruction()
{;}

G4VPhysicalVolume* ExN04DetectorConstruction::Construct()
{
  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
  {
    ExN04Field* myField = new ExN04Field;
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();
    fieldMgr->SetDetectorField(myField);
    fieldMgr->CreateChordFinder(myField);
    fieldIsInitialized = true;
  }

  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

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
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material * Silicon = new G4Material(name="Silicon", z=14., a, density);

  density = 1.032*g/cm3;
  G4Material* Scinti = new G4Material(name="Scintillator", density, nel=2);
  Scinti->AddElement(elC, 9);
  Scinti->AddElement(elH, 10);

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------

  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);
  G4VisAttributes* experimentalHallVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  experimentalHallVisAtt->SetForceWireframe(true);
  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);

  //------------------------------ tracker
  G4VSolid * tracker_tubs
    = new G4Tubs("trkTubs_tubs",trkTubs_rmin,trkTubs_rmax,trkTubs_dz,
                 trkTubs_sphi,trkTubs_dphi);
  G4LogicalVolume * tracker_log
    = new G4LogicalVolume(tracker_tubs,Ar,"trackerT_L",0,0,0);
  G4VPhysicalVolume * tracker_phys
    = new G4PVPlacement(0,G4ThreeVector(),"tracker_phys",tracker_log,
			experimentalHall_phys,false,0);
  G4VisAttributes* tracker_logVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  tracker_logVisAtt->SetForceWireframe(true);
  tracker_log->SetVisAttributes(tracker_logVisAtt);

  //------------------------------ tracker layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * trackerLayer_tubs
    = new G4Tubs("trackerLayer_tubs",trkTubs_rmin,trkTubs_rmax,trkTubs_dz,
                 trkTubs_sphi,trkTubs_dphi);
  G4LogicalVolume * trackerLayer_log
    = new G4LogicalVolume(trackerLayer_tubs,Silicon,"trackerB_L",0,0,0);
  G4VPVParameterisation * trackerParam
    = new ExN04TrackerParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  G4VPhysicalVolume *trackerLayer_phys
    = new G4PVParameterised("trackerLayer_phys",trackerLayer_log,tracker_phys,
			   kXAxis, notrkLayers, trackerParam);
  G4VisAttributes* trackerLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerLayer_logVisAtt->SetForceWireframe(true);
  trackerLayer_log->SetVisAttributes(trackerLayer_logVisAtt);

  //------------------------------ calorimeter
  G4VSolid * calorimeter_tubs
    = new G4Tubs("calorimeter_tubs",caloTubs_rmin,caloTubs_rmax,
		  caloTubs_dz,caloTubs_sphi,caloTubs_dphi);
  G4LogicalVolume * calorimeter_log
    = new G4LogicalVolume(calorimeter_tubs,Scinti,"caloT_L",0,0,0);
  G4VPhysicalVolume * calorimeter_phys
    = new G4PVPlacement(0,G4ThreeVector(),"caloM_P",calorimeter_log,
			experimentalHall_phys,false,0);
  G4VisAttributes* calorimeter_logVisATT
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calorimeter_logVisATT->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeter_logVisATT);

  //------------------------------- Lead layers
  // As an example for Parameterised volume 
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid * caloLayer_tubs
    = new G4Tubs("caloLayer_tubs",caloRing_rmin,caloRing_rmax,
		  caloRing_dz,caloRing_sphi,caloRing_dphi);
  G4LogicalVolume * caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs,Lead,"caloR_L",0,0,0);
  G4VPVParameterisation * calorimeterParam
    = new ExN04CalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  G4VPhysicalVolume * caloLayer_phys
    = new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_phys,
			   kXAxis, nocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_logVisAtt->SetForceWireframe(true);
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------ muon counters
  // As an example of CSG volumes with rotation
  G4VSolid * muoncounter_box
    = new G4Box("muoncounter_box",muBox_width,muBox_thick,
		muBox_length);
  G4LogicalVolume * muoncounter_log
    = new G4LogicalVolume(muoncounter_box,Scinti,"mucounter_L",0,0,0);
  G4VPhysicalVolume * muoncounter_phys;
  for(int i=0; i<nomucounter ; i++)
  {
    G4double phi, x, y, z;
    phi = 360.*deg/nomucounter*i;
    x = muBox_radius*sin(phi);
    y = muBox_radius*cos(phi);
    z = 0.*cm;
    G4RotationMatrix rm;
    rm.rotateZ(phi);
    muoncounter_phys
      = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(x,y,z)),
                          "muoncounter_P",muoncounter_log,
                          experimentalHall_phys,false,i);
  }
  G4VisAttributes* muoncounter_logVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  muoncounter_logVisAtt->SetForceWireframe(true);
  muoncounter_log->SetVisAttributes(muoncounter_logVisAtt);

  //------------------------------------------------------------------
  // Sensitive Detector
  //------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerSDname = "/mydet/tracker";
  ExN04TrackerSD * trackerSD = new ExN04TrackerSD(trackerSDname);
  SDman->AddNewDetector(trackerSD);
  trackerLayer_log->SetSensitiveDetector(trackerSD);

  G4String calorimeterSDname = "/mydet/calorimeter";
  ExN04CalorimeterSD * calorimeterSD = new ExN04CalorimeterSD(calorimeterSDname);
  G4String ROgeometryName = "CalorimeterROGeom";
  G4VReadOutGeometry* calRO = new ExN04CalorimeterROGeometry(ROgeometryName);
  calRO->BuildROGeometry();
  calorimeterSD->SetROgeometry(calRO);
  SDman->AddNewDetector(calorimeterSD);
  calorimeter_log->SetSensitiveDetector(calorimeterSD);

  G4String muonSDname = "/mydet/muon";
  ExN04MuonSD * muonSD = new ExN04MuonSD(muonSDname);
  SDman->AddNewDetector(muonSD);
  muoncounter_log->SetSensitiveDetector(muonSD);

  //------------------------------------------------------------------
  // Digitizer modules
  //------------------------------------------------------------------

  return experimentalHall_phys;
}

