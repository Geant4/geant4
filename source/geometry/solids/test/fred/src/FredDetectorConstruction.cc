//
// FredDetectorConstruction.cc
//
// Implementation of fred's detector
//

#include "FredDetectorConstruction.hh"
#include "FredSensitive.hh"
#include "FredSensMother.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4BREPSolidPCone.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"

#include "../../Boolean/include/G4SubtractionSolid.hh"

#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "globals.hh"

#include <math.h>

//
// Constructor
//
FredDetectorConstruction::FredDetectorConstruction( FredMessenger *ourMessenger ) 
{
	messenger = ourMessenger;
}

FredDetectorConstruction::~FredDetectorConstruction()
{;}


//
// Construct: Build the detector
//
G4VPhysicalVolume* FredDetectorConstruction::Construct()
{
  G4VisAttributes *redStuff = new G4VisAttributes( G4Color(1,0,0) );
	
  //
  // Vaccuum would do fine for materials right now
  //
  G4Material *Vaccuum = new G4Material( "Vaccuum", 18.0, 39.95*g/mole, 1.782e-03*g/cm3 );


  //
  // At the moment, I want something really simple:
  // how about a "hall" containing a single box
  //
	
  /* 
     MEDERNACH Emmanuel
     Aug 2000,

     BEWARE :
     If World is too small for test then Segmentation Fault !!  
  */

  G4Box		  *hallBox = new G4Box( "hall_box", 3*m, 3*m, 3*m );
  G4LogicalVolume	  *hallLog = new G4LogicalVolume( hallBox, Vaccuum, "hall_log", 0, 0, 0 );
  G4VPhysicalVolume *hall	   = new G4PVPlacement( 0, G4ThreeVector(), hallLog, 
							"hall", 0, false, 0 );
						       
  //
  // We usually don't care much about the main volume: just
  // make it invisible
  //
  hallLog->SetVisAttributes( G4VisAttributes::Invisible );
	
  //
  // For the test volume, we have some run-time choices:
  //
  G4RotationMatrix	*rot = new G4RotationMatrix();
	
  G4double startPhi = messenger->StartPhi()*deg,
    deltaPhi = messenger->DeltaPhi()*deg;
  G4int	 numSide  = messenger->NumSide();

  switch( messenger->SelectedVolume() ) {
  case NATALIA: {
    G4double	z_values[3] = { -60.76*mm, -49.14*mm, 102.68*mm };
    G4double	rmin[3]	    = {  6.24*mm, 0*mm, 0*mm },
      rmax[3]     = {  6.24*mm, 6.24*mm, 6.24*mm };
      // Rib thickness 0.41, height 6.42
      startPhi = -atan2( 0.5*0.41, 6.42 );
      deltaPhi = -2.0*startPhi;
      testVolume = new G4Polyhedra( "natalia",
				    startPhi, deltaPhi, 1, 3, z_values, rmin, rmax );
  }
  break;
		
  case VOXEL:
    testVolume = new G4Box( "voxel_test", 1*m, 1*m, 1*m );
    break;
	
  case CONE2:
    /* try to do a cone with a 'pick' */
	  
    testVolume = new G4Cons( "test_cone",
			     1*m, 1.2*m, 0.0*m, 0.2*m, 1*m, startPhi, deltaPhi );
    break;
		
  case CONE:
    testVolume = new G4Cons( "test_cone",
			     1*m, 1.2*m, 0.4*m, 0.6*m, 1*m, startPhi, deltaPhi );
    break;
		
  case TUBS:
    testVolume = new G4Tubs( "test_tube", 1.0*m, 1.2*m, 1*m, startPhi, deltaPhi );
    break;

  case BOX:
    testVolume = new G4Box( "testbox", 1*m, 1*m, 1*m );
    break;
		
  case BOOL1: {
    G4Box	*outside = new G4Box( "testboolout", 1*m, 1*m, 1*m );
    G4Tubs 	*inside = new G4Tubs( "testboolin", 0.0, 0.4*m, 1*m, 0, 360*deg );
    G4Transform3D tran = G4Translate3D( 0.4*m, 0.0, 0.0 );
			
    testVolume = new G4SubtractionSolid( "testbool", (G4VSolid *)outside, (G4VSolid *)inside, tran );
  }
  break;

  case PCON: {
    G4double      z_values[2] = { -1.0*m, 1.0*m };
    G4double      rmin[2]	  = { 1.0*m, 1.2*m },
      rmax[2]	  = { 1.2*m, 1.4*m };
      testVolume = new G4Polycone( "testpcone", 
				   startPhi, deltaPhi, 2, z_values, rmin, rmax );
  }
  break;

  case NEW: {
    /* New volume */
    /*
      You could have:
      G4Sphere.cc  
      G4Trap.cc  
      G4Tubs.cc
      G4CSGSolid.cc  
      G4Para.cc  
      G4Torus.cc   
      G4Trd.cc
    */
    /* OK try a Sphere ..*/
    testVolume = new G4Sphere("testSphere",
			      0.8*m, 1.0*m,
			      startPhi, deltaPhi,
			      0, 2*M_PI);
  }
  break;

  case PCON2: {
    G4double	z_values[5] = { -1.0*m, 0.0*m, 0.0*m, 0.8*m, 1.0*m };
    G4double	rmin[5]	    = {  0.5*m, 0.4*m, 0.0*m, 0.0*m, 0.9*m },
      rmax[5]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polycone( "testpcone", 
				   startPhi, deltaPhi, 5, z_values, rmin, rmax );
  }
  break;
		
  case PCON3: {
    G4double	z_values[8] = { -1.0*m, -0.5*m, -0.5*m, -1.0*m, -1.0*m,  0.7*m,  0.7*m,  1.0*m };
    G4double	rmin[8]     = {  0.6*m,  0.6*m,  0.5*m,  0.5*m,  0.4*m,  0.4*m,  0.4*m,  0.0*m },
      rmax[8]     = {  0.7*m,  0.7*m,  0.8*m,  0.9*m,  1.0*m,  1.0*m,  0.5*m,  0.5*m };
      testVolume = new G4Polycone( "testpcone", 
				   startPhi, deltaPhi, 8, z_values, rmin, rmax );
  }
  break;
		
  case PCON4: {
    double RMINVec[8];
    RMINVec[0] = 30*cm;
    RMINVec[1] = 30*cm;
    RMINVec[2] =  0*cm;
    RMINVec[3] =  0*cm;
    RMINVec[4] =  0*cm; 
    RMINVec[5] =  0*cm;
    RMINVec[6] = 40*cm;
    RMINVec[7] = 40*cm;  

    double RMAXVec[8];
    RMAXVec[0] = 70*cm;
    RMAXVec[1] = 70*cm;
    RMAXVec[2] = 70*cm;
    RMAXVec[3] = 40*cm;
    RMAXVec[4] = 40*cm;
    RMAXVec[5] = 80*cm;
    RMAXVec[6] = 80*cm;
    RMAXVec[7] = 60*cm; 

    double Z_Values[8];
    Z_Values[0] =-20*cm;
    Z_Values[1] =-10*cm;
    Z_Values[2] =-10*cm;
    Z_Values[3] =  0*cm;
    Z_Values[4] = 10*cm;
    Z_Values[5] = 20*cm;
    Z_Values[6] = 30*cm;
    Z_Values[7] = 40*cm;

    testVolume = new G4Polycone ("MyPCone",
				 startPhi       ,
				 deltaPhi     ,
				 8        ,
				 Z_Values ,
				 RMINVec  ,
				 RMAXVec   );
  }
  break;

  case PGON2:
    rot->rotateZ( 360*deg/16 );
  case PGON: {
    G4double	z_values[2] = { -1.0*m, 1.0*m };
    G4double	rmin[2]	    = { 0.8*m, 1.0*m },
      rmax[2]     = { 1.0*m, 1.2*m };
      testVolume = new G4Polyhedra( "testpgon",
				    startPhi, deltaPhi, numSide, 2, z_values, rmin, rmax );
  }
  break;
		
  case PGON3: {
    G4double	z_values[5] = { -1.0*m, 0.0*m, 0.0*m, 0.8*m, 1.0*m };
    G4double	rmin[5]	    = {  0.5*m, 0.4*m, 0.0*m, 0.0*m, 0.9*m },
      rmax[5]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polyhedra( "testpgon",
				    startPhi, deltaPhi, numSide, 5, z_values, rmin, rmax );
  }
  break;
  case PGON4: {
    G4double	z_values[6] = { -0.6*m, 0.0*m,-1.0*m, 0.5*m, 0.5*m, 1.0*m };
    G4double	rmin[6]	    = {  0.5*m, 0.5*m, 0.4*m, 0.4*m, 0.8*m, 0.8*m },
      rmax[6]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polyhedra( "testpgon",
				    startPhi, deltaPhi, numSide, 6, z_values, rmin, rmax );
  }
  break;
  }

  G4LogicalVolume	  *testLog  = new G4LogicalVolume( testVolume, Vaccuum, "test_log", 0, 0, 0 );
  G4VPhysicalVolume *test     = new G4PVPlacement( rot, G4ThreeVector(), testLog, 
						   "test", hallLog, false, 0 );

  //
  // Put some stuff in it, if we want
  //
  if (messenger->SelectedVolume() == VOXEL) {
    G4RotationMatrix	*noRot = new G4RotationMatrix();

    G4Box 			*vxBox = new G4Box( "voxel_x", 0.3*mm, 0.6*m, 0.6*m );
    G4LogicalVolume   	*vxLog1 = new G4LogicalVolume( vxBox, Vaccuum, "x1", 0, 0, 0 );
    G4VPhysicalVolume	*vx1 = new G4PVPlacement( noRot, G4ThreeVector( -0.6*m, 0.0*m, 0.0*m ),
						  vxLog1, "testx1", testLog, false, 0 );
    G4LogicalVolume   	*vxLog2 = new G4LogicalVolume( vxBox, Vaccuum, "x2", 0, 0, 0 );
    G4VPhysicalVolume	*vx2 = new G4PVPlacement( noRot, G4ThreeVector( -0.2*m, 0.0*m, 0.0*m ),
						  vxLog2, "testx2", testLog, false, 0 );
    G4LogicalVolume   	*vxLog3 = new G4LogicalVolume( vxBox, Vaccuum, "x3", 0, 0, 0 );
    G4VPhysicalVolume	*vx3 = new G4PVPlacement( noRot, G4ThreeVector( +0.2*m, 0.0*m, 0.0*m ),
						  vxLog3, "testx3", testLog, false, 0 );
    G4LogicalVolume   	*vxLog4 = new G4LogicalVolume( vxBox, Vaccuum, "x4", 0, 0, 0 );
    G4VPhysicalVolume	*vx4 = new G4PVPlacement( noRot, G4ThreeVector( +0.6*m, 0.0*m, 0.0*m ),
						  vxLog4, "testx4", testLog, false, 0 );

    G4Box 			*vyBox = new G4Box( "voxel_y", 0.8*m, 0.3*mm, 0.6*m );
    G4LogicalVolume   	*vyLog1 = new G4LogicalVolume( vyBox, Vaccuum, "y1", 0, 0, 0 );
    G4VPhysicalVolume	*vy1 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, -0.8*m, 0.0*m ),
						  vyLog1, "testy1", testLog, false, 0 );
    G4LogicalVolume   	*vyLog2 = new G4LogicalVolume( vyBox, Vaccuum, "y2", 0, 0, 0 );
    G4VPhysicalVolume	*vy2 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, -0.7*m, 0.0*m ),
						  vyLog2, "testy2", testLog, false, 0 );
    G4LogicalVolume   	*vyLog3 = new G4LogicalVolume( vyBox, Vaccuum, "y3", 0, 0, 0 );
    G4VPhysicalVolume	*vy3 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, +0.7*m, 0.0*m ),
						  vyLog3, "testy3", testLog, false, 0 );
    G4LogicalVolume   	*vyLog4 = new G4LogicalVolume( vyBox, Vaccuum, "y4", 0, 0, 0 );
    G4VPhysicalVolume	*vy4 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, +0.8*m, 0.0*m ),
						  vyLog4, "testy4", testLog, false, 0 );

    G4Box 			*vzBox = new G4Box( "voxel_z", 0.8*m, 0.8*m, 0.3*mm );
    G4LogicalVolume   	*vzLog1 = new G4LogicalVolume( vzBox, Vaccuum, "z1", 0, 0, 0 );
    G4VPhysicalVolume	*vz1 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, -0.8*m ),
						  vzLog1, "testz1", testLog, false, 0 );
    G4LogicalVolume   	*vzLog2 = new G4LogicalVolume( vzBox, Vaccuum, "z2", 0, 0, 0 );
    G4VPhysicalVolume	*vz2 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, -0.7*m ),
						  vzLog2, "testz2", testLog, false, 0 );
    G4LogicalVolume   	*vzLog3 = new G4LogicalVolume( vzBox, Vaccuum, "z3", 0, 0, 0 );
    G4VPhysicalVolume	*vz3 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, +0.7*m ),
						  vzLog3, "testz3", testLog, false, 0 );
    G4LogicalVolume   	*vzLog4 = new G4LogicalVolume( vzBox, Vaccuum, "z4", 0, 0, 0 );
    G4VPhysicalVolume	*vz4 = new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, +0.8*m ),
						  vzLog4, "testz4", testLog, false, 0 );
  }

  //
  // Red seems an appropriate color
  //
  testLog->SetVisAttributes( redStuff );

  //
  // Too simple?? Yeah. Let's make our test volume sensitive.
	
  G4SDManager *sensitiveMan = G4SDManager::GetSDMpointer();
	
  FredSensitive *sensitive = new FredSensitive( "/fred/test" );
  sensitiveMan->AddNewDetector( sensitive );
  testLog->SetSensitiveDetector( sensitive );
	
  //
  // And while we're at it, do the same to mother volume
	
  FredSensMother *sensMother = new FredSensMother( "/fred/mother" );
  sensitiveMan->AddNewDetector( sensMother );
  hallLog->SetSensitiveDetector( sensMother );
	
  //
  // Tell our "messenger" about this test volume
  //
  messenger->SetTestVolume( testVolume );
	
  return hall;
}
