// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "ExN06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"

ExN06DetectorConstruction::ExN06DetectorConstruction()
{
  expHall_x = 10.*m;
  expHall_y = 10.*m;
  expHall_z = 10.*m;

  tank_x = 5.*m;
  tank_y = 5.*m;
  tank_z = 5.*m;

  bubble_x = 0.5*m;
  bubble_y = 0.5*m;
  bubble_z = 0.5*m;
}

ExN06DetectorConstruction::~ExN06DetectorConstruction(){;}

G4VPhysicalVolume* ExN06DetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4double a, z, density;
  G4String name, symbol;
  G4int nel;

// Air
// ---

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

// Water
// -----

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);

  density = 1.0*g/cm3;
  G4Material* Water = new G4Material(name="Water", density, nel=2);

  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

/////////////////////////////////////////////
// Generate & Add Material Properties Table
/////////////////////////////////////////////

  const G4int NUMENTRIES = 32;

  G4double PPCKOV[NUMENTRIES] =
            { 2.034E-9*GeV, 2.068E-9*GeV, 2.103E-9*GeV, 2.139E-9*GeV,
              2.177E-9*GeV, 2.216E-9*GeV, 2.256E-9*GeV, 2.298E-9*GeV,
              2.341E-9*GeV, 2.386E-9*GeV, 2.433E-9*GeV, 2.481E-9*GeV,
              2.532E-9*GeV, 2.585E-9*GeV, 2.640E-9*GeV, 2.697E-9*GeV,
              2.757E-9*GeV, 2.820E-9*GeV, 2.885E-9*GeV, 2.954E-9*GeV,
              3.026E-9*GeV, 3.102E-9*GeV, 3.181E-9*GeV, 3.265E-9*GeV,
              3.353E-9*GeV, 3.446E-9*GeV, 3.545E-9*GeV, 3.649E-9*GeV,
              3.760E-9*GeV, 3.877E-9*GeV, 4.002E-9*GeV, 4.136E-9*GeV };

//            { 2.038E-9*GeV, 2.072E-9*GeV, 2.107E-9*GeV, 2.143E-9*GeV,
//              2.181E-9*GeV, 2.220E-9*GeV, 2.260E-9*GeV, 2.302E-9*GeV,
//              2.346E-9*GeV, 2.391E-9*GeV, 2.438E-9*GeV, 2.486E-9*GeV,
//              2.537E-9*GeV, 2.590E-9*GeV, 2.645E-9*GeV, 2.702E-9*GeV,
//              2.763E-9*GeV, 2.825E-9*GeV, 2.891E-9*GeV, 2.960E-9*GeV,
//              3.032E-9*GeV, 3.108E-9*GeV, 3.188E-9*GeV, 3.271E-9*GeV,
//              3.360E-9*GeV, 3.453E-9*GeV, 3.552E-9*GeV, 3.656E-9*GeV,
//              3.767E-9*GeV, 3.884E-9*GeV, 4.010E-9*GeV, 4.144E-9*GeV };

  G4double RINDEX1[NUMENTRIES] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};

//            { 1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33,
//              1.33, 1.33, 1.34, 1.34, 1.34, 1.34, 1.34,
//              1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
//              1.34, 1.34, 1.35, 1.35, 1.35, 1.35, 1.35,
//              1.35, 1.35, 1.35, 1.35 };

  G4double RINDEX2[NUMENTRIES] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4double ABSORPTION1[NUMENTRIES] =
           {344.8*cm,  408.2*cm,  632.9*cm,  917.4*cm, 1234.6*cm, 1388.9*cm,
           1515.2*cm, 1724.1*cm, 1886.8*cm, 2000.0*cm, 2631.6*cm, 3571.4*cm,
           4545.5*cm, 4761.9*cm, 5263.2*cm, 5263.2*cm, 5555.6*cm, 5263.2*cm,
           5263.2*cm, 4761.9*cm, 4545.5*cm, 4166.7*cm, 3703.7*cm, 3333.3*cm,
           3000.0*cm, 2850.0*cm, 2700.0*cm, 2450.0*cm, 2200.0*cm, 1950.0*cm,
           1750.0*cm, 1450.0*cm };

  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", PPCKOV, RINDEX1, NUMENTRIES);
  myMPT1->AddProperty("ABSLENGTH",PPCKOV, ABSORPTION1, NUMENTRIES);
  Water->SetMaterialPropertiesTable(myMPT1);

  G4MaterialPropertiesTable *myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PPCKOV, RINDEX2, NUMENTRIES);
  Air->SetMaterialPropertiesTable(myMPT2);

//	------------- Volumes --------------

//	The experimental Hall
//	---------------------

  G4Box * expHall_box
    = new G4Box("World",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume * expHall_log
    = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);

  G4VPhysicalVolume * expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"World",expHall_log,0,false,0);

//	The Water Tank
//	--------------

  G4Box * waterTank_box
    = new G4Box("Tank",tank_x,tank_y,tank_z);

  G4LogicalVolume * waterTank_log
    = new G4LogicalVolume(waterTank_box,Water,"Tank",0,0,0);

  G4RotationMatrix *rot1=new G4RotationMatrix();
  rot1->rotateZ(M_PI*0.125);

//  G4VPhysicalVolume * waterTank_phys
//    = new G4PVPlacement(rot1,G4ThreeVector(),"Tank",
//			waterTank_log,expHall_phys,false,0);   

  G4VPhysicalVolume * waterTank_phys
    = new G4PVPlacement(0,G4ThreeVector(),"Tank",
                        waterTank_log,expHall_phys,false,0);

//      The Air Bubble
//      --------------

  G4Box * bubbleAir_box
    = new G4Box("Bubble",bubble_x,bubble_y,bubble_z);

  G4LogicalVolume * bubbleAir_log
    = new G4LogicalVolume(bubbleAir_box,Air,"Bubble",0,0,0);

  G4RotationMatrix *rot2=new G4RotationMatrix();
  rot2->rotateZ(M_PI*0.25);

//  G4VPhysicalVolume * bubbleAir_phys
//    = new G4PVPlacement(rot2,G4ThreeVector(0,2.5*m,0),"Bubble",
//                        bubbleAir_log,waterTank_phys,false,0);

  G4VPhysicalVolume * bubbleAir_phys
    = new G4PVPlacement(0,G4ThreeVector(0,2.5*m,0),"Bubble",
                        bubbleAir_log,waterTank_phys,false,0);

//	------------- Surfaces --------------

  G4OpticalSurface * OpWaterSurface =
                                new G4OpticalSurface("WaterSurface");

  G4LogicalBorderSurface * WaterSurface = 
	new G4LogicalBorderSurface("WaterSurface",
				    waterTank_phys,expHall_phys,
				    OpWaterSurface);

  OpWaterSurface->SetType(dielectric_metal);
  OpWaterSurface->SetFinish(polished);
  OpWaterSurface->SetModel(glisur);

  if( WaterSurface->GetVolume1() == waterTank_phys ) G4cout << " Equal " << endl;
  if( WaterSurface->GetVolume2() == expHall_phys   ) G4cout << " Equal " << endl;

  G4OpticalSurface * OpAirSurface =
                              new G4OpticalSurface("AirSurface");

  G4LogicalSkinSurface * AirSurface = 
	new G4LogicalSkinSurface("AirSurface",
				  bubbleAir_log,
				  OpAirSurface);

  OpAirSurface->SetType(dielectric_dielectric);
  OpAirSurface->SetFinish(ground);
  OpAirSurface->SetModel(unified);

  if( AirSurface->GetLogicalVolume() == bubbleAir_log ) G4cout << " Equal " << endl;

  G4LogicalBorderSurface * Tmp1Surface = WaterSurface->
				 GetSurface(waterTank_phys,expHall_phys);

//  if (Tmp1Surface == *WaterSurface) G4cout << " Equal "     << endl;

  G4LogicalSkinSurface * Tmp2Surface = AirSurface->GetSurface(bubbleAir_log);

//  if (Tmp2Surface == *AirSurface  ) G4cout << " Equal "     << endl;

  G4OpticalSurface * TmpOpSurface = Tmp2Surface->GetOpticalSurface();
  TmpOpSurface->DumpInfo();

/////////////////////////////////////////////
// Generate & Add Material Properties Table
/////////////////////////////////////////////

  const G4int NUM = 2;

  G4double PP[NUM] =
            { 2.038E-9*GeV, 4.144E-9*GeV };

  G4double RINDEX[NUM] =
            { 1.35, 1.40 };
  G4double SPECULARLOBECONSTANT[NUM] =
            { 0.3, 0.3 };
  G4double SPECULARSPIKECONSTANT[NUM] =
            { 0.2, 0.2 };
  G4double BACKSCATTERCONSTANT[NUM] =
            { 0.2, 0.2 };

  G4MaterialPropertiesTable *myST1 = new G4MaterialPropertiesTable();

  myST1->AddProperty("RINDEX", PP, RINDEX, NUM);
  myST1->
  AddProperty("SPECULARLOBECONSTANT", PP, SPECULARLOBECONSTANT, NUM);
  myST1->
  AddProperty("SPECULARSPIKECONSTANT", PP, SPECULARSPIKECONSTANT, NUM);
  myST1->
  AddProperty("BACKSCATTERCONSTANT", PP, BACKSCATTERCONSTANT, NUM);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);

  G4double REFLECTIVITY[NUM] =
            { 0.3, 0.5 };
  G4double EFFICIENCY[NUM] =
            { 0.8, 1.0 };

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", PP, REFLECTIVITY, NUM);
  myST2->AddProperty("EFFICIENCY", PP, EFFICIENCY, NUM);

  OpAirSurface->SetMaterialPropertiesTable(myST2);

  return expHall_phys;
}
