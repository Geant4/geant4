// $Id: Tst10DetectorConstruction.cc,v 1.2 1999-12-15 14:54:42 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "Tst10DetectorConstruction.hh"
#include "Tst10DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst10DetectorConstruction::Tst10DetectorConstruction() {
  detectorMessenger = new Tst10DetectorMessenger (this);
}

Tst10DetectorConstruction::~Tst10DetectorConstruction(){
;
}

void Tst10DetectorConstruction::SwitchDetector( void ) {
  G4UImanager::GetUIpointer()->ApplyCommand("/run/geometryModified"); 
}
G4VPhysicalVolume* Tst10DetectorConstruction::SelectDetector( G4String val ) {

//------------------- A Volume ----------------------

if (val == "Sphere")
// 	aVolume = new G4Sphere ( "aSphere", 8.0*cm, 10.0*cm, 
// 					 0.0*deg, 300.0*deg,30.0*deg, 110.0*deg);
	aVolume = new G4Sphere ( "aSphere", 8.0*cm, 10.0*cm, 
					 0.0*deg, 360.0*deg,0.0*deg, 130.0*deg);
else if (val == "Box") 				 
	aVolume = new G4Box ( "aBox", 10*cm, 10*cm, 10*cm );
else if (val == "Cone") 			 
	aVolume = new G4Cons ( "aCone", 2*cm, 6*cm, 8*cm, 14*cm, 10*cm, 10*deg,
          	300*deg );	
else if (val == "Tube")
	aVolume = new G4Tubs ( "aTube", 5*cm, 10*cm, 7*cm, 70*deg, 100*deg);
else if (val == "Hype")
	aVolume = new G4Hype ("aHype", 10*cm, 20*cm, 0*deg, 360*deg, 10*cm );
else if (val == "Torus")
	aVolume = new G4Torus ("aTorus", 10*cm, 15*cm, 20*cm, 0*deg, 60*deg);
else if (val == "Para")
	aVolume = new G4Para ("aPara", 8*cm, 10*cm, 12*cm, 30*deg, 45*deg, 60*deg);
else if (val == "Trd")
	aVolume = new G4Trd ("aTrd", 8*cm, 10*cm, 7*cm, 9*cm, 10*cm);
else {
  cout << "You don't select valid shape " << G4endl;
	exit (1);
}

  G4Box * Hall 
       = new G4Box("Hall", 1*m, 1*m, 1*m );
  G4LogicalVolume * Hall_log 
      = new G4LogicalVolume (Hall, Water, "Hall_L", 0,0,0);
  PhysicalVolume
      = new G4PVPlacement(0,G4ThreeVector(),"Hall_P", 
	                  Hall_log, 0, false, 0);
   
  G4LogicalVolume * aVolume_log 
		= new G4LogicalVolume(aVolume, Water1, "aVolume_L", 0,0,0);
  G4VPhysicalVolume * aVolume_phys1
		= new G4PVPlacement(0,G4ThreeVector(50*cm, 0*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys2
		= new G4PVPlacement(0,G4ThreeVector(-50*cm, 0*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys3
		= new G4PVPlacement(0,G4ThreeVector(0*cm, 50*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys4
		= new G4PVPlacement(0,G4ThreeVector(0*cm, -50*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys5
		= new G4PVPlacement(0,G4ThreeVector(0*cm, 0*cm, 50*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys6
		= new G4PVPlacement(0,G4ThreeVector(0*cm, 0*cm, -50*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);


// ------------ Surfaces definition ------------------

	G4LogicalBorderSurface* BorderSurfaces[12]; 	
	BorderSurfaces[0] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys1,
                               aSurface);
	BorderSurfaces[1] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys2,
                               aSurface);
	BorderSurfaces[2] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys3,
                               aSurface);
	BorderSurfaces[3] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys4,
                               aSurface);
	BorderSurfaces[4] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys5,
                               aSurface);
	BorderSurfaces[5] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys6,
                               aSurface);
	BorderSurfaces[6] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys1,
                               PhysicalVolume,
                               aSurface);
	BorderSurfaces[7] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys2,
                               PhysicalVolume,
                               aSurface);
	BorderSurfaces[8] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys3,
                               PhysicalVolume,
                               aSurface);
	BorderSurfaces[9] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys4,
                               PhysicalVolume,
                               aSurface);
	BorderSurfaces[10] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys5,
                               PhysicalVolume,
                               aSurface);
	BorderSurfaces[11] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys6,
                               PhysicalVolume,
                               aSurface);
															 
G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()
                ->SetWorldVolume(PhysicalVolume);
//if (G4VVisManager::GetConcreteInstance()!=NULL)
//  G4VVisManager::GetConcreteInstance()->SetWorldVolume(PhysicalVolume);

cout << "You select " << val << " detector" << G4endl;

return PhysicalVolume;
}

void Tst10DetectorConstruction::SetMaterial( void ) {

  G4String name, symbol;
  G4double density = 1.00*g/cm3;
  G4double a, iz;
  Water = new G4Material(name="Water", density, 2);
  Water1 = new G4Material(name="Water1", density, 2);
  a = 1*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  Water->AddElement(elH, .66);
  Water->AddElement(elO, .34);
  Water1->AddElement(elH, .66);
  Water1->AddElement(elO, .34);

  const G4int NUMENTRIES = 5;
  G4double RINDEX_WATER [NUMENTRIES];
  G4double RINDEX_WATER1 [NUMENTRIES];
  G4double REFLECTIVITY [NUMENTRIES];
  G4double EFFICIENCY [NUMENTRIES];
	
  for (int i=0; i<NUMENTRIES; i++) {
	  RINDEX_WATER1[i]=5.0;
	  RINDEX_WATER[i]=1.33;
	  REFLECTIVITY[i]=0.9;
	  EFFICIENCY[i]=1.0;
	}	
  G4double PHENERGY[NUMENTRIES] =
            { 0.0, 1.0, 2.0, 3.0, 4.0};
  G4MaterialPropertiesTable *WaterMPT = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *WaterMPT1 = new G4MaterialPropertiesTable();
  WaterMPT->AddProperty("RINDEX", PHENERGY, RINDEX_WATER, NUMENTRIES);
  WaterMPT1->AddProperty("RINDEX", PHENERGY, RINDEX_WATER1, NUMENTRIES);
  Water->SetMaterialPropertiesTable(WaterMPT);
  Water1->SetMaterialPropertiesTable(WaterMPT1);

  aSurface = new G4OpticalSurface ( "aSurface" );
  aSurface->SetType(dielectric_metal);
  aSurface->SetFinish(polishedfrontpainted);
  aSurface->SetModel(glisur);  
  G4MaterialPropertiesTable* SurfaceMPT = new G4MaterialPropertiesTable();
  SurfaceMPT->AddProperty("REFLECTIVITY", PHENERGY, REFLECTIVITY, NUMENTRIES);
  SurfaceMPT->AddProperty("EFFICIENCY", PHENERGY, EFFICIENCY, NUMENTRIES);
  aSurface->SetMaterialPropertiesTable ( SurfaceMPT );


}
G4VPhysicalVolume* Tst10DetectorConstruction::Construct() {

  SetMaterial();

//-------------------Hall ----------------------------------
	
  G4VPhysicalVolume * Hall_phys = SelectDetector ("Sphere");

  return Hall_phys;
}


