// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		Tst18GeometryConstruction.cc
//
// Version:		0.b.3
// Date:		30/06/99
// Author:		P R Truscott
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 30 June 1999, P R Truscott, DERA UK
// Version number update 0.b.2 -> 0.b.3, but no functional change.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "Tst18GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "PhysicalConstants.h"
////////////////////////////////////////////////////////////////////////////////
//
Tst18GeometryConstruction::Tst18GeometryConstruction ()
{
  //
  //
  // Define size of experiment hall and volumes which will fill it.
  //
  universe_x = 500.0*km;
  universe_y = 500.0*km;
  universe_z = 500.0*km;

  aSphere_r1 = 250.0*km;
  aSphere_r2 = 125.0*km;
}
////////////////////////////////////////////////////////////////////////////////
//
Tst18GeometryConstruction::~Tst18GeometryConstruction ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* Tst18GeometryConstruction::Construct ()
{
  //
  //
  // Define materials.
  //
  G4double a, iz, z, density,pressure, temperature;
  G4String name, symbol;

 
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* Vacuum   = new G4Material("Vacuum",
                                        1., 1.01*g/mole, density,
                                        kStateGas,temperature,pressure);    
  a = 26.98154*g/mole;
  density = 2.70*g/cm3;
  G4Material* Aluminium = new G4Material(name="aluminium", z=13., a, density);

  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element(name="silicon", symbol="Si", iz=14., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 2.65*g/cm3;
  G4Material* SiliconDioxide =
    new G4Material(name="silicon oxide", density, nel=2);
  SiliconDioxide->AddElement(elSi, 1);
  SiliconDioxide->AddElement(elO,  2);
  //
  //
  // Define bodies, logical volumes and physical volumes.
  // First define the experimental hall.
  //
  G4Box * universe_box
    = new G4Box("universe_b",universe_x,universe_y,universe_z);
  G4LogicalVolume * universe_log
    = new G4LogicalVolume(universe_box,Vacuum,"universe_L",0,0,0);
  G4VPhysicalVolume * universe_phys
    = new G4PVPlacement(0,G4ThreeVector(),"universe_P",
                        universe_log,0,false,0);

  //
  //
  // Define outer sphere.
  //
  G4Sphere * aSphere_sph1
    = new G4Sphere("aSphere1",0,aSphere_r1,0.*deg,360.*deg,0.*deg,180.*deg);
  G4LogicalVolume * aSphere_log1
    = new G4LogicalVolume(aSphere_sph1,Aluminium,"aSphere_L1",0,0,0);
  G4VPhysicalVolume * aSphere_phys1
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"aSphere_P1",
    aSphere_log1,universe_phys,false,0);

  //
  //
  // Define middle sphere.
  //
  G4Sphere * aSphere_sph2
    = new G4Sphere("aSphere2",0,aSphere_r2,0.*deg,360.*deg,0.*deg,180.*deg);
  G4LogicalVolume * aSphere_log2
    = new G4LogicalVolume(aSphere_sph2,SiliconDioxide,"aSphere_L2",0,0,0);
  G4VPhysicalVolume * aSphere_phys2
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"aSphere_P2",aSphere_log2,
                      universe_phys,false,0);

  return universe_phys;
}

