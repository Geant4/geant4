// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildShapes.cc,v 1.2 1999-12-15 14:54:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "BuildShapes.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"

#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"

// #define USE_BREP_SOLID

#if defined USE_BREP_SOLID

#include"G4BREPSolidPCone.h"
#include"G4BREPSolidPolyhedra.h"

#endif

G4VPhysicalVolume* TheWorld ()
{
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  G4Box *myWorldBox= new G4Box("WBox",10.0, 10.0, 10.0);  // Big enough!
  G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,Air,
                                                    "WLog", 0, 0, 0);
  myWorldLog -> SetVisAttributes (&G4VisAttributes::Invisible);
  G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
                                                 "WPhys",
                                                 myWorldLog,
                                                 0,false,0);
  return myWorldPhys;
}

G4VPhysicalVolume* BuildBox()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,0,0));

		//----- Solid
	G4Box *box = new G4Box ("Box", 1 ,4, 9 );

		//----- Logical Volume
	G4LogicalVolume *boxLog=new G4LogicalVolume( box, &myMaterial,
						    "BoxLog",0,0,0);
	boxLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *boxPhys=new G4PVPlacement( 0,G4ThreeVector(),
						  "BoxPhys",
						  boxLog,
						  theWorld,
						  false,0);

	return theWorld;

} // BuildBox


G4VPhysicalVolume* BuildCylinder()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(0,1,0));

		//----- Solid
	G4Tubs *tubs = new G4Tubs ("Tubs", 0. , 1.5, 2.0 , 0., 2.*M_PI);

		//----- Logical Volume
	G4LogicalVolume *tubsLog=new G4LogicalVolume( tubs, &myMaterial,
						    "tubsLog",0,0,0);
	tubsLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *tubsPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "tubsPhys",
						   tubsLog,
						   theWorld,
						   false,0);

	return theWorld;

}// BuildCyllinder



G4VPhysicalVolume* BuildTubs()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(0,1,0));

		//----- Solid
	G4Tubs *tubs = new G4Tubs ("Tubs", 1.0 , 1.5, 2.0 , 1.5, 4.5 );

		//----- Logical Volume
	G4LogicalVolume *tubsLog=new G4LogicalVolume( tubs, &myMaterial,
						    "tubsLog",0,0,0);
	tubsLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *tubsPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "tubsPhys",
						   tubsLog,
						   theWorld,
						   false,0);

	return theWorld;

}// BuildTubs



G4VPhysicalVolume* BuildCons()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(0,0,1));

		//----- Solid
	G4Cons *cons = new G4Cons ("Cons", 1.0 , 1.5, 0.5, 1.0 ,2.0, 1.5, 4.5 );

		//----- Logical Volume
	G4LogicalVolume *consLog=new G4LogicalVolume( cons, &myMaterial,
						    "consLog",0,0,0);
	consLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *consPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "consPhys",
						   consLog,
						   theWorld,
						   false,0);

	return theWorld;

}// BuildCons


G4VPhysicalVolume* BuildTrd()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(0,1,1));

		//----- Solid
	G4Trd *trd = new G4Trd ("Trd", 1.0 , 0.5, 1.0, 0.5, 4.0 );

		//----- Logical Volume
	G4LogicalVolume *trdLog=new G4LogicalVolume( trd, &myMaterial,
						    "trdLog",0,0,0);
	trdLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *trdPhys=new G4PVPlacement( 0,G4ThreeVector(),
						  "trdPhys",
						  trdLog,
						  theWorld,
						  false,0);

	return theWorld;

}// BuildTrd



G4VPhysicalVolume* BuildTrap()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,0,1));

		//----- Solid
	G4Trap *trap = new G4Trap ("Trap", 2.0, 0.26, 0.0, 1.0, \
				   0.7, 1.2, 0.0, 0.5, 0.35, 0.6, 0.0 );

		//----- Logical Volume
	G4LogicalVolume *trapLog=new G4LogicalVolume( trap, &myMaterial,
						    "trapLog",0,0,0);
	trapLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *trapPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "trapPhys",
						   trapLog,
						   theWorld,
						   false,0);

	return theWorld;

}// BuildTrap


G4VPhysicalVolume* BuildSphereFull()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,1,0));

		//----- Solid
	G4Sphere *Sphere = new G4Sphere ("Sphere", 1.6, 1.8, 0.0, (2.0*M_PI), 0.0, (2.0*M_PI) );

		//----- Logical Volume
	G4LogicalVolume *SphereLog=new G4LogicalVolume( Sphere, &myMaterial,
						    "SphereLog",0,0,0);
	SphereLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *SpherePhys=new G4PVPlacement( 0,G4ThreeVector(),
						     "SpherePhys",
						     SphereLog,
						     theWorld,
						     false,0);

	return theWorld;

} // BuildSphereFull


G4VPhysicalVolume* BuildSphereSeg()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,1,0));

		//----- Solid
	G4Sphere *Sphere = new G4Sphere ("Sphere", 1.6, 1.8, 0.0, 6.0, 0.5, 2.5 );

		//----- Logical Volume
	G4LogicalVolume *SphereLog=new G4LogicalVolume( Sphere, &myMaterial,
						    "SphereLog",0,0,0);
	SphereLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *SpherePhys=new G4PVPlacement( 0,G4ThreeVector(),
						     "SpherePhys",
						     SphereLog,
						     theWorld,
						     false,0);

	return theWorld;

} // BuildSphereSeg



G4VPhysicalVolume* BuildPara()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,0,0));

		//----- Solid
//	G4Para *para = new G4Para ("Para", 1.0 ,1.0, 3.0, 0.0, 0.0, 0.0 );
	G4Para *para = new G4Para ("Para", 1.0 ,1.0, 3.0, M_PI/6.0, M_PI/12.0, M_PI/12.0 );

		//----- Logical Volume
	G4LogicalVolume *paraLog=new G4LogicalVolume( para, &myMaterial,
						    "ParaLog",0,0,0);
	paraLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *paraPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "ParaPhys",
						   paraLog,
						   theWorld,
						   false,0);

	return theWorld;

} // BuildPara




G4VPhysicalVolume* BuildPCon()
{
#if defined USE_BREP_SOLID
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- shape
	G4double sphi = 1.8 ;
	G4double dphi = 4.71 ;
	G4int    nz   = 4    ;
	G4double DzArray   [] = { -3.0, 1.0, 1.0, 3.0 };
	G4double RminArray [] = {  1.0, 0.5, 1.5, 1.0 };
	G4double RmaxArray [] = {  1.5, 1.0, 2.0, 1.5 };

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(0,1,0));

		//----- Solid
	G4BREPSolidPCone *pcon = new G4BREPSolidPCone ("Pcon", sphi , dphi, nz, DzArray[0], \
						     DzArray, RminArray, RmaxArray     );

		//----- Logical Volume
	G4LogicalVolume *pconLog=new G4LogicalVolume( pcon, &myMaterial,
						    "PconLog",0,0,0);
	pconLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *pconPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "PconPhys",
						   pconLog,
						   theWorld,
						   false,0);

	return theWorld;

#else 
	return NULL ;
#endif

} // BuildPCon


G4VPhysicalVolume* BuildPGon()
{
#if defined USE_BREP_SOLID
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- shape
	G4double sphi = 1.8 ;
	G4double dphi = 4.71 ;
	G4int    nside   = 5    ;
	G4int    nz   = 4    ;
	G4double DzArray   [] = { -3.0, 1.0, 1.0, 3.0 };
	G4double RminArray [] = {  1.0, 0.5, 1.5, 1.0 };
	G4double RmaxArray [] = {  1.5, 1.0, 2.0, 1.5 };

		//----- Vis Attributes
	G4VisAttributes* pVA  = new G4VisAttributes (G4Colour(1,0,0));

		//----- Solid
	G4BREPSolidPolyhedra *pgon = new G4BREPSolidPolyhedra ("Pcon", sphi , dphi, nside, nz, \
	                                                       DzArray[0], \
						               DzArray, RminArray, RmaxArray     );

		//----- Logical Volume
	G4LogicalVolume *pgonLog=new G4LogicalVolume( pgon, &myMaterial,
						    "PgonLog",0,0,0);
	pgonLog -> SetVisAttributes ( pVA );

		//----- physical volume
	G4PVPlacement *pgonPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "PgonPhys",
						   pgonLog,
						   theWorld,
						   false,0);
	
	return theWorld;
#else
	return NULL ;
#endif
} // BuildPGon

G4VPhysicalVolume* BuildForcedWireframeBox()
{
        G4VPhysicalVolume* theWorld = TheWorld ();

		//----- material
	G4double z = 13.;
	G4double a = 26.98*g/mole;
	G4double density = 2.7*g/cm3;
	G4Material myMaterial("Aluminium", z, a, density);
	
		//----- Vis Attributes
	G4VisAttributes* pVA   = new G4VisAttributes (G4Colour(1,0,0)); // for box
	pVA->SetForceWireframe( true );
	G4VisAttributes* pVA2  = new G4VisAttributes (G4Colour(0,1,1)); // for tubs

		//----- Solid
	G4Box *box   = new G4Box  ("Box", 2 ,4, 6 );
	G4Tubs *tubs = new G4Tubs ("Tubs", 1.0 , 1.5, 2.0 , 1.5, 4.5 );

		//----- Logical Volume
	G4LogicalVolume *boxLog=new G4LogicalVolume( box, &myMaterial,
						    "BoxLog",0,0,0);
	boxLog -> SetVisAttributes ( pVA );
	G4LogicalVolume *tubsLog=new G4LogicalVolume( tubs, &myMaterial,
						    "tubsLog",0,0,0);
	tubsLog -> SetVisAttributes ( pVA2 );

		//----- physical volume
	G4PVPlacement *boxPhys=new G4PVPlacement( 0,G4ThreeVector(),
						  "BoxPhys",
						  boxLog,
						  theWorld,
						  false,0);
	G4PVPlacement *tubsPhys=new G4PVPlacement( 0,G4ThreeVector(),
						   "tubsPhys",
						   tubsLog,
						   boxPhys ,false,0);
	return theWorld;

} // G4VPhysicalVolume* BuidlForcedWireFrameBox()
