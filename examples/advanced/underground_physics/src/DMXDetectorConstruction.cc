//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
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
#include "G4SphericalSurface.hh"

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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()  
{
  theUserLimits  = NULL; 
  fUseUserLimits = false;
  //  theMaxTimeCuts = 1000000. * s;
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts = DBL_MAX;
  theMaxStepSize = DBL_MAX;

  theRoomTimeCut = 10000. * s;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction() {;}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials() {

  G4double density,      // density
    a,                   // atomic mass
    z;                   // atomic number
  G4String name,         // name
    symbol;              // symbol
  G4int ncomponents,     // n components
    iz,                  // number of protons
    in;                  // number of nuceons
  G4double abundance,    // abundance
    temperature,         // temperature
    pressure;            // pressure


  // making vacuum
  G4Material* vacuum = new G4Material 
    (name="Vacuum", z=1., a=1.*g/mole, density=1.e-20*g/cm3,
     kStateGas, temperature=0.1*kelvin, pressure=1.e-20*bar);


  // xenons
  G4Element* elementXe = new G4Element( "Xenon", "Xe", 54., 131.29*g/mole );
  G4Material* LXe = new G4Material
     ("LXe", 3.02*g/cm3, 1, kStateLiquid, 173.15*kelvin, 1.5*atmosphere );
  G4Material* GXe = new G4Material
     ("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere );
  LXe->AddElement( elementXe, 1);
  GXe->AddElement( elementXe, 1);

  const G4int NUMENTRIES = 2;
  G4double LXe_PP[NUMENTRIES]    = { 7.07*eV, 7.07*eV };
  G4double LXe_SCINT[NUMENTRIES] = { 1.0, 1.0 };
  G4double LXe_RIND[NUMENTRIES]  = { 1.57, 1.57 };
  G4double LXe_ABSL[NUMENTRIES]  = { 35.*cm, 35.*cm};
  G4MaterialPropertiesTable *LXe_mt = new G4MaterialPropertiesTable();
  LXe_mt->AddProperty("SCINTILLATION", LXe_PP, LXe_SCINT, NUMENTRIES);
  LXe_mt->AddProperty("RINDEX",        LXe_PP, LXe_RIND,  NUMENTRIES);
  LXe_mt->AddProperty("ABSLENGTH",     LXe_PP, LXe_ABSL,  NUMENTRIES);
  LXe->SetMaterialPropertiesTable(LXe_mt);

  G4double GXe_PP[NUMENTRIES]    = { 7.07*eV, 7.07*eV };
  G4double GXe_SCINT[NUMENTRIES] = { 1.0, 1.0 };
  G4double GXe_RIND[NUMENTRIES]  = { 1.00, 1.00 };
  G4double GXe_ABSL[NUMENTRIES]  = { 100*m, 100*m};
  G4MaterialPropertiesTable *GXe_mt = new G4MaterialPropertiesTable();
  GXe_mt->AddProperty("SCINTILLATION", GXe_PP, GXe_SCINT, NUMENTRIES);
  GXe_mt->AddProperty("RINDEX",        GXe_PP, GXe_RIND,  NUMENTRIES);
  GXe_mt->AddProperty("ABSLENGTH",     GXe_PP, GXe_ABSL,  NUMENTRIES);
  GXe->SetMaterialPropertiesTable(GXe_mt);


  // making quartz
  G4Element* O  = new G4Element
    (name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
  G4Element* Si = new G4Element
    (name="Silicon",symbol="Si" , z= 14., a=28.09*g/mole);
  G4Material* quartz = new G4Material
    (name="quartz", density=2.200*g/cm3, ncomponents=2);
  quartz->AddElement(Si, 1);
  quartz->AddElement(O , 2);

  G4double quartz_PP[NUMENTRIES]   = { 6.69*eV, 7.50*eV }; // lambda range 4 ri
  G4double quartz_RIND[NUMENTRIES] = { 1.575, 1.628 };     // ref index
  G4double quartz_ABSL[NUMENTRIES] = { 1.5*cm, 1.5*cm };   // atten length
  G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();
  quartz_mt->AddProperty("RINDEX", quartz_PP, quartz_RIND, NUMENTRIES);
  quartz_mt->AddProperty("ABSLENGTH", quartz_PP, quartz_ABSL, NUMENTRIES);
  quartz->SetMaterialPropertiesTable(quartz_mt);


  // aluminium
  G4Element* Al = new G4Element
    (name="Aluminium"  ,symbol="Al" , z= 13., a=26.98*g/mole);  
  G4Material* metalAl = new G4Material
    (name="MetalAluminium", density=2.700*g/cm3, ncomponents=1);
  metalAl->AddElement(Al, 1);


  // iron
  G4Element* Fe = new G4Element
    (name="Iron"  ,symbol="Fe" , z= 26., a=55.85*g/mole);  
  G4Material* metalFe = new G4Material
    (name="MetalIron", density=7.874*g/cm3, ncomponents=1);
  metalFe->AddElement(Fe, 1);


  // stainless steel
  G4Element* C  = new G4Element( "Carbon", "C",   6. , 12.011*g/mole);
  G4Element* Co = new G4Element( "Cobalt", "Co", 27. , 58.9332*g/mole);
  G4Material* ssteel = new G4Material
    (name="Steel", density=7.7*g/cm3, ncomponents=3);
  ssteel->AddElement(C, 0.04);
  ssteel->AddElement(Fe, 0.88);
  ssteel->AddElement(Co, 0.08);


  // copper
  G4Element* Cu = new G4Element
    (name="Copper"  ,symbol="Cu" , z= 29., a=63.55*g/mole);  
  G4Material* metalCu = new G4Material
    (name="MetalCopper", density=8.960*g/cm3, ncomponents=1);
  metalCu->AddElement(Cu, 1);

  // lead
  G4Element* Pb = new G4Element
    (name="Lead",symbol="Pb" , z= 82., a=207.2*g/mole);
  G4Material* metalPb = new G4Material
    (name="MetalLead", density=11.340*g/cm3, ncomponents=1);
  metalPb->AddElement(Pb, 1);


  // Americium:
  G4Isotope* Am241 = new G4Isotope
    (name="Americium241", iz= 95, in=241, a=241.0*g/mole);
  G4Element* Am = new G4Element
    (name="Americium241", "Am", ncomponents=1);
  Am->AddIsotope(Am241, abundance=1);
  G4Material* sourceAm = new G4Material
    (name="AmericiumSource", density=13.61*g/cm3, ncomponents=1);
  sourceAm->AddElement(Am, 1);

  // air
  G4Element* N = new G4Element
    (name="Nitrogen",symbol="N" , z= 7., a=14.00674*g/mole);
  G4Material* Air = new G4Material
    ("AIR", 1.2929*kg/m3, 2, kStateGas, 300.00*kelvin, 1.0*atmosphere);
  Air->AddElement(N, 0.8);
  Air->AddElement(O , 0.2);


  //concrete
  G4Element* H = new G4Element
    (name="Hydrogen",symbol="H" , z= 1., a=1.00794*g/mole);
  G4Element* Ca = new G4Element
    (name="Calcium",symbol="Ca" , z= 20., a=40.078*g/mole);
  G4Material* concrete = new G4Material
    (name="Concrete", density=2.3*g/cm3, ncomponents=6);
  concrete->AddElement(Si, 0.227915);
  concrete->AddElement(O, 0.60541);
  concrete->AddElement(H, 0.09972);
  concrete->AddElement(Ca, 0.04986);
  concrete->AddElement(Al, 0.014245);
  concrete->AddElement(Fe, 0.00285);

  //water
  G4Material* water = new G4Material
    (name="water", density=1.00*g/cm3, ncomponents=2);
  water->AddElement(H , 2);
  water->AddElement(O , 1);

  // print materials
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  //  G4cout << *(G4Isotope::GetIsotopeTable())   << G4endl;
  //  G4cout << *(G4Element::GetElementTable())   << G4endl;

  // assign materials
       world_mat = concrete;
         lab_mat = Air;
      jacket_mat = ssteel;
      vacuum_mat = vacuum;
      vessel_mat = ssteel;
    detector_mat = GXe;
    CuShield_mat = metalCu;
    liqPhase_mat = LXe;
       alpha_mat = metalPb;
      recess_mat = metalPb;
   americium_mat = sourceAm;
        ring_mat = ssteel;
      mirror_mat = metalAl;
        grid_mat = LXe;
         pmt_mat = quartz;
      phcath_mat = metalAl;
}


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineField() {

  G4double EField = 0.0*kilovolt/cm;

  // create electric field
  G4UniformElectricField* elecField = 
     new G4UniformElectricField(G4ThreeVector(0,EField,0));

  // equation of motion
  G4EqMagElectricField* EquationOfMotion = 
    new G4EqMagElectricField(elecField);

  // stepper for equation of motion
  G4MagIntegratorStepper* DMXStepper = 
    new G4ClassicalRK4(EquationOfMotion); 

  // chordfinder
  G4ChordFinder* DMXChordFinder =
    new G4ChordFinder(elecField, 1.0e-3*mm, DMXStepper);

  // field manager
  G4FieldManager* DMXFieldManager = new G4FieldManager();
  DMXFieldManager->SetChordFinder(DMXChordFinder);  
  G4TransportationManager::GetTransportationManager()
    -> SetFieldManager(DMXFieldManager);
}
*/


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.75, .75, .75) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;
  //  G4Colour  green   (0.0, 1.0, 0.0) ;
  //  G4Colour  lgreen  (0.0, .75, 0.0) ;



  // Universe - room wall - CONCRETE ************************************

  G4double worldWidth  = 5.0*m;
  G4double worldLength = 6.0*m;
  G4double worldHeight = 3.0*m;

  G4Box* world_box= new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log = new G4LogicalVolume(world_box, world_mat, "world_log");
  world_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "world_phys", world_log, NULL, false,0);

  G4VisAttributes* world_vat= new G4VisAttributes(white);
  //  world_log->SetVisAttributes(G4VisAttributes::Invisible);
  world_vat->SetVisibility(true);
  world_log->SetVisAttributes(world_vat);


  // Lab Space - AIR *****************************************************

  G4double wallThick = 30.*cm;
  G4double labWidth  = worldWidth  - wallThick;
  G4double labLength = worldLength - wallThick;
  G4double labHeight = worldHeight - wallThick;

  G4Box* lab_box= new G4Box
     ("lab_box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight );
  lab_log = new G4LogicalVolume(lab_box, lab_mat, "lab_log");
  lab_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lab_phys", 
     lab_log, world_phys, false,0);

  G4VisAttributes* lab_vat= new G4VisAttributes(white);
  lab_log->SetVisAttributes(G4VisAttributes::Invisible);
  //  lab_log->SetVisibility(true);
  lab_log->SetVisAttributes(lab_vat);


  // outer vacuum jacket volume: stainless steel *************************

  G4double jacketRadius     = 15.0*cm;
  G4double jacketHeight     = 40.0*cm;
  G4double jacketMetalThick = 3.0*mm;

  G4Tubs* jacket_tube=new G4Tubs("jacket_tube",
     0.*cm, jacketRadius, 0.5*jacketHeight, 0.*deg, 360.*deg);
  jacket_log = new G4LogicalVolume(jacket_tube, jacket_mat, "jacket_log");
  jacket_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "jacket_phys", jacket_log, lab_phys, false,0);

  G4VisAttributes* jacket_vat= new G4VisAttributes(grey);
  //  jacket_log->SetVisAttributes(G4VisAttributes::Invisible);
  jacket_log->SetVisAttributes(jacket_vat);


  // vacuum **************************************************************

  G4double vacuumRadius = jacketRadius - jacketMetalThick;
  G4double vacuumHeight = jacketHeight - jacketMetalThick;

  G4Tubs* vacuum_tube=new G4Tubs("vacuum_tube",
     0.*cm, vacuumRadius, 0.5*vacuumHeight, 0.*deg, 360.*deg);
  vacuum_log = new G4LogicalVolume(vacuum_tube, vacuum_mat, "vacuum_log");
  vacuum_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "vacuum_phys", vacuum_log, jacket_phys, false,0);

  //  G4VisAttributes* vacuum_vat= new G4VisAttributes(lgrey);
  vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);


  // inner vessel jacket volume: stainless steel ************************

  G4double vesselRadius = 10.0*cm;
  G4double vesselHeight = 30.0*cm;
  G4double vesselMetalThick = 3.0*mm;

  G4Tubs* vessel_tube=new G4Tubs("vessel_tube",
     0.*cm, vesselRadius, 0.5*vesselHeight, 0.*deg, 360.*deg);
  vessel_log = new G4LogicalVolume(vessel_tube, vessel_mat, "vessel_log");
  vessel_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), 
     "vessel_phys", vessel_log, vacuum_phys, false,0);

  G4VisAttributes* vessel_vat= new G4VisAttributes(grey);
  //  vessel_log->SetVisAttributes(G4VisAttributes::Invisible);
  vessel_log->SetVisAttributes(vessel_vat);


  // *********************************************************************
  // grid#1 to mirror surface: 21.75 mm
  // LXe height = 15.75 mm, gXe height = 6.00 mm
  // NB: Increased liquid height by 1mm - to take away problem with 
  // over-lapping volumes/ring pronounced from liquid phase..........
  // *********************************************************************

  // detector volume: gas phase ******************************************

  G4double fullDetectorRadius = vesselRadius - vesselMetalThick;
  G4double fullDetectorHeight = vesselHeight - vesselMetalThick;

  G4Tubs* detector_tube=new G4Tubs("detector_tube",
     0.*cm, fullDetectorRadius, 0.5*fullDetectorHeight, 0.*deg, 360.*deg);
  detector_log = new G4LogicalVolume
     (detector_tube, detector_mat, "detector_log");
  detector_phys=new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), 
     "detector_phys", detector_log, vessel_phys, false,0);

  G4VisAttributes* detector_vat= new G4VisAttributes(yellow);
  //  detector_log->SetVisAttributes(G4VisAttributes::Invisible);
  detector_log->SetVisAttributes(detector_vat);


  // liquid phase *******************************************************

  G4double detectorRadius     = 4.0*cm;
  G4double detectorHeight     = 21.3*cm;
  G4double liqPhaseHeight    = 20.40*cm;
  G4double liqPhaseRadius    = fullDetectorRadius;
  G4double liqPhaseVPosition = -0.5*(detectorHeight-liqPhaseHeight);

  G4Tubs* liqPhase_tube=new G4Tubs("liqPhase_tube", 0.*cm, 
     liqPhaseRadius, 0.5*liqPhaseHeight, 0.*deg, 360.*deg);
  liqPhase_log =
     new G4LogicalVolume(liqPhase_tube, liqPhase_mat, "liqPhase_log");
  liqPhase_phys=new G4PVPlacement(0, 
    G4ThreeVector(0.*cm, 0.*cm, liqPhaseVPosition), 
    "liqPhase_phys", liqPhase_log, detector_phys, false, 0);

  // attributes
  G4VisAttributes* liqPhase_vat= new G4VisAttributes(yellow);
  liqPhase_vat->SetVisibility(true);
  liqPhase_log->SetVisAttributes(liqPhase_vat);


  // Cu Shield **********************************************************

  G4double CuShieldHeight      = 17.7*cm;
  G4double CuShieldThickness   = 2.4*mm;
  G4double CuShieldOuterRadius = 3.0*cm;
  G4double CuShieldInnerRadius = CuShieldOuterRadius-CuShieldThickness;
  G4double CuShieldVPosition   = -0.5*(liqPhaseHeight-CuShieldHeight);

  G4Tubs* CuShield_tube=new G4Tubs("CuShield_tube", CuShieldInnerRadius,
     CuShieldOuterRadius, 0.5*CuShieldHeight, 0.*deg, 360.*deg);
  CuShield_log = 
     new G4LogicalVolume(CuShield_tube, CuShield_mat, "CuShield_log");
  CuShield_phys=new G4PVPlacement(0, 
     G4ThreeVector(0.*cm, 0.*cm, CuShieldVPosition), 
     "CuShield_phys", CuShield_log, liqPhase_phys, false, 0);
  G4VisAttributes* CuShield_vat= new G4VisAttributes(magenta);
  CuShield_vat->SetVisibility(true);
  CuShield_log->SetVisAttributes(CuShield_vat);

  // Cu shield surface
  G4OpticalSurface* OpCuShieldSurface = new G4OpticalSurface
    ("ShieldSurface", unified, polished, dielectric_metal);
  G4LogicalBorderSurface* ShieldSurface;
  ShieldSurface = new G4LogicalBorderSurface
    ("Shield", liqPhase_phys, CuShield_phys, OpCuShieldSurface);

  const G4int NUM = 2;
  G4double CuShield_PP[NUM]   = { 7.0*eV, 7.50*eV };
  G4double CuShield_REFL[NUM] = { 0.271, 0.230 };
  G4MaterialPropertiesTable *CuShield_mt = new G4MaterialPropertiesTable();
  CuShield_mt->AddProperty("REFLECTIVITY", CuShield_PP, CuShield_REFL, NUM);
  OpCuShieldSurface->SetMaterialPropertiesTable(CuShield_mt);


  // rings ***************************************************************

  G4double ringHeight      =  4.*mm;
  G4double ringOuterRadius =  detectorRadius;
  G4double ringInnerRadius =  CuShieldOuterRadius;
  G4double ringVOffset     =  0.5*ringHeight;
  G4double ringVPosition   =  0.5*detectorHeight-ringVOffset;

  G4Tubs* ring_tube=new G4Tubs("ring_tube", ringInnerRadius,
     ringOuterRadius, 0.5*ringHeight, 0.*deg, 360.*deg);
  ring_log = new G4LogicalVolume(ring_tube, ring_mat, "ring_log");

  // optical surface: ring materials table
  G4double ring_PP[NUM]   = { 7.00*eV, 7.50*eV };
  G4double ring_REFL[NUM] = { 0.25, 0.25 };
  G4MaterialPropertiesTable *ring_mt = new G4MaterialPropertiesTable();
  ring_mt->AddProperty("REFLECTIVITY", ring_PP, ring_REFL, NUM);

  G4OpticalSurface* OpRingSurface = new G4OpticalSurface("RingSurface", 
     unified, polished, dielectric_metal);
  // omitted last argument which is surface roughness if it's non-polished 
  // - i.e. ground
  OpRingSurface->SetMaterialPropertiesTable(ring_mt);

  // rings
  ring_phys_gas[0]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition),
     "ring_phys0",ring_log,detector_phys,false, 0);
  G4LogicalBorderSurface* RingSurface_gas0;
  RingSurface_gas0 = new G4LogicalBorderSurface
    ("Ring", detector_phys, ring_phys_gas[0], OpRingSurface);

  ring_phys_gas[1]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.0*mm),
     "ring_phys1",ring_log, detector_phys, false, 0);
  G4LogicalBorderSurface* RingSurface_gas1;
  RingSurface_gas1 = new G4LogicalBorderSurface
    ("Ring", detector_phys, ring_phys_gas[1], OpRingSurface);


  // LIQUID Phase starts here:
  ringVPosition+=0.5*(detectorHeight-liqPhaseHeight);

  ring_phys_liq[0]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight),
     "ring_phys2",ring_log,liqPhase_phys, false, 0);
  G4LogicalBorderSurface* RingSurface_liq0;
  RingSurface_liq0 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[0], OpRingSurface);

  ring_phys_liq[1]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys3",ring_log, liqPhase_phys, false, 0);
  G4LogicalBorderSurface* RingSurface_liq1;
  RingSurface_liq1 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[1], OpRingSurface);

  ring_phys_liq[2]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight),
     "ring_phys4",ring_log, liqPhase_phys, false, 0);
  G4LogicalBorderSurface* RingSurface_liq2;
  RingSurface_liq2 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[2], OpRingSurface);

  ring_phys_liq[3]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight),
     "ring_phys5",ring_log, liqPhase_phys, false, 0);
  G4LogicalBorderSurface* RingSurface_liq3;
  RingSurface_liq3 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[3], OpRingSurface);

  ring_phys_liq[4]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys6",ring_log, liqPhase_phys,false, 0);
  G4LogicalBorderSurface* RingSurface_liq4;
  RingSurface_liq4 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[4], OpRingSurface);

  ring_phys_liq[5]=new G4PVPlacement(0,
     G4ThreeVector(0.*cm, 0.*cm, ringVPosition-=ringHeight+1.75*mm),
     "ring_phys7",ring_log, liqPhase_phys,false, 0);
  G4LogicalBorderSurface* RingSurface_liq5;
  RingSurface_liq5 = new G4LogicalBorderSurface
    ("Ring", liqPhase_phys, ring_phys_liq[5], OpRingSurface);


  G4VisAttributes* ring_vat= new G4VisAttributes(lgrey);
  ring_vat->SetVisibility(true);
  ring_log->SetVisAttributes(ring_vat);


  // Mirror *************************************************************

  G4double mirrorHeight    = 2.0*mm;
  G4double mirrorRadius    = ringInnerRadius;
  G4double mirrorVOffset   = 0.5*ringHeight;
  G4double mirrorVPosition = 0.5*detectorHeight-mirrorVOffset;

  G4Tubs* mirror_tube=new G4Tubs("mirror_tube", 0.*cm, mirrorRadius,
     0.5*mirrorHeight, 0.*deg, 360.*deg);
  mirror_log =
    new G4LogicalVolume(mirror_tube, mirror_mat, "mirror_log");
  mirror_phys=new G4PVPlacement(0, 
     G4ThreeVector(0.*cm, 0.*cm, mirrorVPosition),
     "mirror_phys", mirror_log, detector_phys, false, 0);

  G4VisAttributes* mirror_vat= new G4VisAttributes(red);
  mirror_vat->SetVisibility(true);
  //  mirror_vat->SetForceSolid(true);
  mirror_log->SetVisAttributes(mirror_vat);


  // mirror surface
  G4OpticalSurface * OpMirrorSurface = new G4OpticalSurface
    ("MirrorSurface", unified, polished, dielectric_metal);
  G4LogicalBorderSurface* MirrorSurface;
  MirrorSurface = new G4LogicalBorderSurface
    ("Mirror", detector_phys, mirror_phys, OpMirrorSurface);

  G4double mirror_PP[NUM]   = { 7.00*eV, 7.50*eV };
  G4double mirror_REFL[NUM] = { 0.70, 0.70 };
  G4MaterialPropertiesTable *mirror_mt = new G4MaterialPropertiesTable();
  mirror_mt->AddProperty("REFLECTIVITY", mirror_PP, mirror_REFL, NUM);
  OpMirrorSurface->SetMaterialPropertiesTable(mirror_mt);


  // Grids  *************************************************************

  G4double gridHeight     = 0.100*mm;
  G4double gridRadius     = ringInnerRadius;
  G4double grid1VOffset   = 5.5*ringHeight+2.75*mm;
  G4double grid1VPosition = 0.5*liqPhaseHeight+(detectorHeight-liqPhaseHeight)
                            - grid1VOffset;
  G4double grid2VOffset   = 6.5*ringHeight+4.50*mm;
  G4double grid2VPosition = 0.5*liqPhaseHeight+(detectorHeight-liqPhaseHeight)
                            - grid2VOffset;

  G4Tubs* grid_tube=new G4Tubs("grid_tube", 0.*cm, gridRadius,
     0.5*gridHeight, 0.*deg, 360.*deg);

  grid1_log = new G4LogicalVolume(grid_tube, grid_mat, "grid1_log");
  grid1_phys=new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, grid1VPosition),
     "grid1_phys", grid1_log, liqPhase_phys, false, 0);
  grid2_log = new G4LogicalVolume(grid_tube, grid_mat, "grid2_log");
  grid2_phys=new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, grid2VPosition),
     "grid2_phys", grid2_log, liqPhase_phys, false, 0);

  G4VisAttributes* grid_vat= new G4VisAttributes(red);
  grid_vat->SetVisibility(true);
  grid1_log->SetVisAttributes(grid_vat);
  grid2_log->SetVisAttributes(grid_vat);
  
  
  // alpha source holder ************************************************
  
  G4double alphaHeight     = 1.0*mm; // position where alpha starts
  G4double recessHeight    = 0.23*mm;  // totals lead thickness = 1.23 mm
  G4double alphaRadius     = 0.65*mm;  
  G4double recessRadius    = 0.40*mm;

  G4double alphaVOffset    = grid1VOffset-0.5*alphaHeight;
  G4double alphaVPosition  = 0.5*liqPhaseHeight+(detectorHeight-liqPhaseHeight)
                             - alphaVOffset;
  G4double recessVOffset   = 0.5*(alphaHeight - recessHeight);
  G4double recessVPosition = alphaVPosition + recessVOffset;

  G4Tubs* alpha_tube  = new G4Tubs("alpha_tube", 0.*cm, alphaRadius,
     0.5*alphaHeight,  0.*deg, 360.*deg);
  G4Tubs* recess_tube = new G4Tubs("recess_tube", recessRadius, alphaRadius,
     0.5*recessHeight, 0.*deg, 360.*deg);

  alpha_log = new G4LogicalVolume(alpha_tube, alpha_mat, "alpha_log");
  alpha_phys  = new G4PVPlacement(0, G4ThreeVector(0., 0., alphaVPosition),
    "alpha_phys", alpha_log, liqPhase_phys, false, 0);

  recess_log = new G4LogicalVolume(recess_tube, recess_mat, "recess_log");
  recess_phys  = new G4PVPlacement(0, G4ThreeVector(0., 0., recessVPosition),
    "recess_phys", recess_log, liqPhase_phys, false, 0);

  G4VisAttributes* alpha_vat= new G4VisAttributes(white);
  alpha_vat->SetVisibility(true);
  alpha_log ->SetVisAttributes(alpha_vat);
  recess_log ->SetVisAttributes(alpha_vat);

  // alpha source HOLDER surface

  G4OpticalSurface* OpAlphaSurface = new G4OpticalSurface("AlphaSurface", 
     unified, polished, dielectric_metal);
  G4LogicalBorderSurface* AlphaSurface;
  AlphaSurface = new G4LogicalBorderSurface
    ("Alpha", liqPhase_phys, alpha_phys, OpAlphaSurface);

  G4OpticalSurface* OpRecessSurface = new G4OpticalSurface("RecessSurface", 
     unified, polished, dielectric_metal);
  G4LogicalBorderSurface* RecessSurface;
  RecessSurface = new G4LogicalBorderSurface
    ("Recess", liqPhase_phys, recess_phys, OpRecessSurface);

  G4double alpha_PP[NUM]   = { 7.00*eV, 7.50*eV };
  G4double alpha_REFL[NUM] = { 0.1, 0.1 };
  G4MaterialPropertiesTable *alpha_mt = new G4MaterialPropertiesTable();
  alpha_mt->AddProperty("REFLECTIVITY", alpha_PP, alpha_REFL, NUM);
  OpAlphaSurface->SetMaterialPropertiesTable(alpha_mt);

  G4double recess_PP[NUM]   = { 7.00*eV, 7.50*eV };
  G4double recess_REFL[NUM] = { 0.1, 0.1 };
  G4MaterialPropertiesTable *recess_mt = new G4MaterialPropertiesTable();
  recess_mt->AddProperty("REFLECTIVITY", recess_PP, recess_REFL, NUM);
  OpRecessSurface->SetMaterialPropertiesTable(recess_mt);

  // americium ***********************************************************

  G4double americiumHeight    = 20.*nanometer;
  G4double americiumRadius    = recessRadius;
  //  G4double americiumVOffset   = 0.5*(alphaHeight-americiumHeight)-recessHeight;
  G4double americiumVOffset   = 0.5*(alphaHeight-americiumHeight);
  G4double americiumVPosition = americiumVOffset;

  sourceZ = -0.5*(detectorHeight-liqPhaseHeight)+
    alphaVPosition + americiumVPosition;
  G4cout << "Calibration source: Z= " << sourceZ/mm << " mm" << G4endl;

  G4Tubs* americium_tube = new G4Tubs("americium_tube", 0.*cm,
     americiumRadius, 0.5*americiumHeight, 0.*deg, 360.*deg);
  americium_log = new G4LogicalVolume(americium_tube, americium_mat,
     "americium_log");
  americium_phys = new G4PVPlacement(0, G4ThreeVector(0., 0.,
     americiumVPosition),"americium_phys", americium_log, alpha_phys,false,0);

  G4VisAttributes* americium_vat= new G4VisAttributes(cyan);
  americium_vat->SetVisibility(true);
  americium_vat->SetForceSolid(true);
  americium_log->SetVisAttributes(americium_vat);


  // Photomultiplier: ETL 9829 QA ****************************************

  G4double pmtHeight    = 12.0*cm;
  G4double pmtRadius    = 2.6*cm;
  G4double pmtVOffset   = 5.0*cm;
  G4double pmtVPosition = -0.5*(liqPhaseHeight-pmtHeight)+pmtVOffset;

  //  G4double windowVOffset   = 0.5*pmtHeight - 2.*pmtRadius*cos(30.0*deg);
  //  G4double windowVPosition = pmtVPosition + windowVOffset;

  G4Sphere* pmt_window = new G4Sphere("pmt_sphere", 0.*cm, 2.*pmtRadius, 
     0.*deg, 360.*deg, 0.*deg, 30.0*deg);
  G4Tubs* pmt_tube=new G4Tubs("pmt_tube", 0.*cm,  pmtRadius, 0.5*pmtHeight,
     0.*deg, 360.*deg); 

  G4RotationMatrix rotMatrixpmt;       // unit rotation matrix
  G4double anglepmt_Phys = 0.0*deg;    // rotational angle
  rotMatrixpmt.rotateY(anglepmt_Phys); // rot matrix
  
  G4UnionSolid* pmt_sol = new G4UnionSolid("pmt_sol", pmt_tube, pmt_window,
    G4Transform3D(rotMatrixpmt, G4ThreeVector(0,0,0.5*pmtHeight
    -2.*pmtRadius*cos(30.0*deg))));

  pmt_log = new G4LogicalVolume(pmt_sol, pmt_mat, "pmt_log");
  pmt_phys = new G4PVPlacement(0,G4ThreeVector(0.*cm, 0.*cm, pmtVPosition),
     "pmt_phys", pmt_log, liqPhase_phys, false, 0);

  G4OpticalSurface* pmt_opsurf = new G4OpticalSurface("pmt_opsurf",
     unified, polished, dielectric_dielectric);
  G4LogicalBorderSurface* pmt_surf;
  pmt_surf = new G4LogicalBorderSurface
    ("pmt_surf", liqPhase_phys, pmt_phys, pmt_opsurf);

  G4VisAttributes* pmt_vat= new G4VisAttributes(blue);
  pmt_vat->SetForceSolid(true);
  pmt_vat->SetVisibility(true);
  pmt_log->SetVisAttributes(pmt_vat);


  // photocathode *******************************************************

  //  G4double phcathRadius      = 22.5*mm;
  //  G4double phcathVOffset     = 2.*pmtRadius*cos(30.*deg);
  G4double phcathVPosition   = 0.*cm;

  G4Sphere* phcath_sol = new G4Sphere("phcath_sphere",
     2.*pmtRadius-1.6*mm, 2.*pmtRadius-1.59*mm, 0.*deg, 360.*deg, 0.*deg, 
     27.0*deg);

  phcath_log  = new G4LogicalVolume(phcath_sol, phcath_mat, "phcath_log");
  phcath_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., phcathVPosition),
     "phcath_phys", phcath_log, pmt_phys, false, 0);

  G4OpticalSurface*  phcath_opsurf = new G4OpticalSurface("phcath_opsurf",
     unified, polished, dielectric_dielectric);
  G4LogicalBorderSurface* phcath_surf;
  phcath_surf = new G4LogicalBorderSurface
    ("phcath_surf", pmt_phys, phcath_phys, phcath_opsurf);

  G4double phcath_PP[NUM]   = { 7.00*eV, 7.50*eV };
  G4double phcath_REFL[NUM] = { 0.0, 0.0};
  G4MaterialPropertiesTable* phcath_mt = new G4MaterialPropertiesTable();
  phcath_mt->AddProperty("REFLECTIVITY", phcath_PP, phcath_REFL, NUM);
  phcath_opsurf->SetMaterialPropertiesTable(phcath_mt);

  G4VisAttributes* phcath_vat= new G4VisAttributes(lblue);
  phcath_vat->SetForceSolid(true);
  phcath_vat->SetVisibility(true);
  phcath_log->SetVisAttributes(phcath_vat);


  // ......................................................................
  // attach user limits ...................................................

  G4double RoomTimeCut = 10000. * s;

  world_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,RoomTimeCut));
  lab_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,RoomTimeCut));
  jacket_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,RoomTimeCut));
  vacuum_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,RoomTimeCut));
  vessel_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,RoomTimeCut));
  detector_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  liqPhase_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  CuShield_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  ring_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  mirror_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  grid1_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  grid2_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  alpha_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  recess_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  americium_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  pmt_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));
  phcath_log->SetUserLimits
    (new G4UserLimits(theMaxStepSize,DBL_MAX,theMaxTimeCuts));


  // ......................................................................
  // sensitive detectors ..................................................

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String name="/DMXDet/LXeSD";
  LXeSD = new DMXScintSD(name, this);
  SDman->AddNewDetector(LXeSD);
  liqPhase_log->SetSensitiveDetector(LXeSD);

  SDman = G4SDManager::GetSDMpointer();
  name="/DMXDet/pmtSD";
  pmtSD = new DMXPmtSD(name, this);
  SDman->AddNewDetector(pmtSD);
  phcath_log->SetSensitiveDetector(pmtSD);


  return world_phys;

}

