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
// FredDetectorConstruction.cc
//
// Implementation of fred's detector
//

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4Hype.hh"
#include "G4QuadrangularFacet.hh"
#include "G4Tet.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"

#include "G4SubtractionSolid.hh"

#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "globals.hh"

#include <cmath>

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
// Private methods
//
G4ExtrudedSolid*  FredDetectorConstruction::CreateExtrudedSolid1() const
{
  // Extruded solid with triangular polygon
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(  0.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  

  return new G4ExtrudedSolid("test_xtru1", polygon, 30.*cm, 
                             G4TwoVector(), 1.0, G4TwoVector(), 1.0);
}                             

G4ExtrudedSolid*  FredDetectorConstruction::CreateExtrudedSolid2() const
{
  // Box defined as Extruded solid
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  
  return new G4ExtrudedSolid("test_xtru2", polygon, 30.*cm, 
                             G4TwoVector(), 1.0, G4TwoVector(), 1.0);
}

G4ExtrudedSolid*  FredDetectorConstruction::CreateExtrudedSolid3() const
{
  // Extruded solid with 4 z-sections
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));
  
  std::vector<G4ExtrudedSolid::ZSection> zsections;
  zsections.push_back(G4ExtrudedSolid::ZSection(-40.*cm, G4TwoVector(-20.*cm, 10.*cm), 1.5));
  zsections.push_back(G4ExtrudedSolid::ZSection( 10.*cm, G4TwoVector(  0.*cm,  0.*cm), 0.5));
  zsections.push_back(G4ExtrudedSolid::ZSection( 15.*cm, G4TwoVector(  0.*cm,  0.*cm), 0.7));
  zsections.push_back(G4ExtrudedSolid::ZSection( 40.*cm, G4TwoVector( 20.*cm, 20.*cm), 0.9));

  return new G4ExtrudedSolid("test_xtru3", polygon, zsections);
}

G4ExtrudedSolid*  FredDetectorConstruction::CreateExtrudedSolid4() const
{
  // Another extruded solid, where polygon decomposition was failing
  // in Geant4 9.1
  std::vector<G4TwoVector> polygon; 
  polygon.push_back( G4TwoVector(-20.*cm,  10.*cm) );
  polygon.push_back( G4TwoVector(-20.*cm,  25.*cm) );
  polygon.push_back( G4TwoVector( 10.*cm,  25.*cm) );
  polygon.push_back( G4TwoVector( 10.*cm, -10.*cm) );
  polygon.push_back( G4TwoVector( 20.*cm, -10.*cm) );
  polygon.push_back( G4TwoVector( 20.*cm, -25.*cm) );
  polygon.push_back( G4TwoVector(-10.*cm, -25.*cm) );
  polygon.push_back( G4TwoVector(-10.*cm,  10.*cm) );
  
  return new G4ExtrudedSolid("test_xtru3", polygon, 20.*cm, 
                             G4TwoVector(), 1.0, G4TwoVector(), 1.0);
}

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

  G4Box		  *hallBox = new G4Box( "hall_box", 8*m, 8*m, 8*m );
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

  switch ( messenger->SelectedVolume() ) {
  
    //
    // special tests
    //

    case NATALIA: {
      G4double	z_values[3] = { -60.76*mm, -49.14*mm, 102.68*mm };
      G4double	rmin[3]	    = {  6.24*mm, 0*mm, 0*mm };
      G4double	rmax[3]     = {  6.24*mm, 6.24*mm, 6.24*mm };
      // Rib thickness 0.41, height 6.42
      startPhi = -std::atan2( 0.5*0.41, 6.42 );
      deltaPhi = -2.0*startPhi;
      testVolume = new G4Polyhedra("natalia",
                                    startPhi, deltaPhi, 1, 3, z_values, rmin, rmax );
      }                                    
      break;

    case VOXEL:
      testVolume = new G4Box( "test_voxel", 1*m, 1*m, 1*m );
      break;
	
  //
  // CSG solids
  //
  
    case BOX:
      testVolume = new G4Box( "test_box", 1*m, 1*m, 1*m );
      break;
		
    case CONE:
      testVolume = new G4Cons("test_cone",
      			      1*m, 1.2*m, 0.4*m, 0.6*m, 1*m, startPhi, deltaPhi );
      // SBT test - case c
      // testVolume = new G4Cons("test_cone",
      //                              0.0*m, 1.0*m, 0.5*m, 1.0*m, 1*m, 0.0, 360.0*deg );
      // SBT test - case d
      // testVolume = new G4Cons("test_cone",
      //                              0.0*m, 1.0*m, 0.0*m, 1.0*m, 1*m, 0.0, 90.0*deg );
      // SBT test - case e
      // testVolume = new G4Cons("test_cone",
      //                              0.0*m, 1.0*m, 0.0*m, 1.0*m, 1*m, 20.0*deg, 181.0*deg );
      // SBT test - case f
      // testVolume = new G4Cons("test_cone",
      //                              0.5*m, 1.0*m, 0.7*m, 1.2*m, 1*m, 20.0*deg, 350.0*deg );
      // SBT test - case g
      // testVolume = new G4Cons("test_cone",
      //                              0.0*m, 0.2*m, 0.8*m, 1.0*m, 0.0001*m, 10.0*deg, 90.0*deg );
     fprintf(stderr,"OK defining a Cone \n");
      break;

    case CONE2:
      // try to do a cone with a 'pick'
      testVolume = new G4Cons( "test_cone2",
			       1*m, 1.2*m, 0.0*m, 0.2*m, 1*m, startPhi, deltaPhi );
      fprintf(stderr,"OK defining a Cone2 \n");
      break;
		
    case ORB:
      testVolume = new G4Orb ("test_orb", 1.0*m);
      fprintf(stderr,"OK defining an Orb\n");
      break;
    
    case PARA:
      // SBT test case b
      testVolume = new G4Para("test_para",
                                     1.0*m, 1.0*m, 1.0*m, 30.0*deg, 0.0*deg, 0.0*deg);
      // SBT test case c
      // testVolume = new G4Para("test_para",
      //			       1.0*m, 1.0*m, 1.0*m, 30.0*deg, 30.0*deg, 0.0*deg);
      // SBT test case d
      // testVolume = new G4Para("test_para",
      //			       1.0*m, 1.0*m, 1.0*m, 30.0*deg, 30.0*deg, 30.0*deg);
      // SBT test case e
      // testVolume = new G4Para("test_para",
      //			       0.001*m, 1.0*m, 2.0*m, 30.0*deg, 30.0*deg, 30.0*deg);
      fprintf(stderr,"OK defining a Para \n");
      break;

    case SPHERE:
      testVolume = new G4Sphere ("test_sphere", 0.8*m, 1.0*m, startPhi, deltaPhi, 0.0, pi);
      // SBT test case a
      // testVolume = new G4Sphere ("test_sphere", 0.0*m, 1.0*m, 0.0*deg, 360.0*deg, 0.0*deg, 180.0*deg);
      // SBT test case b
      // testVolume = new G4Sphere ("test_sphere", 0.5*m, 1.0*m, 0.0*deg, 360.0*deg, 0.0*deg, 180.0*deg);
      // SBT test case c
      // testVolume = new G4Sphere ("test_sphere", 0.0*m, 1.0*m, 0.0*deg, 90.0*deg, 0.0*deg, 180.0*deg);
      // SBT test case d
      // testVolume = new G4Sphere ("test_sphere", 0.5*m, 1.0*m, 0.0*deg, 90.0*deg, 0.0*deg, 180.0*deg);
      // SBT test case e
      // testVolume = new G4Sphere ("test_sphere", 0.0*m, 1.0*m, 0.0*deg, 360.0*deg, 0.0*deg, 90.0*deg);
      // SBT test case e
      // testVolume = new G4Sphere ("test_sphere", 0.5*m, 1.0*m, 0.0*deg, 360.0*deg, 0.0*deg, 90.0*deg);
      // SBT test case f
      // testVolume = new G4Sphere ("test_sphere", 0.0*m, 1.0*m, 0.0*deg, 90.0*deg, 0.0*deg, 90.0*deg);
      // SBT test case g
      // testVolume = new G4Sphere ("test_sphere", 0.5*m, 1.0*m, 0.0*deg, 90.0*deg, 0.0*deg, 90.0*deg);
      fprintf(stderr,"OK defining a Sphere \n");
      break;
      
    case TORUS1:
      testVolume = new G4Torus("test_torus1",
                                      0.2*m, 0.4*m, 1.2*m, startPhi, deltaPhi);
      // SBT test case a
      // testVolume = new G4Torus("test_torus1",
      //			        0.0*m, 0.4*m, 1.0*m, 0.0*deg, 360.0*deg);
      // SBT test case b
      // testVolume = new G4Torus("test_torus1",
      //			        0.2*m, 0.4*m, 1.0*m, 0.0*deg, 360.0*deg);
      // SBT test case c
      // testVolume = new G4Torus("test_torus1",
      //			        0.0*m, 0.4*m, 1.0*m, 0.0*deg, 90.0*deg);
      // SBT test case d
      // testVolume = new G4Torus("test_torus1",
      //			        0.2*m, 0.4*m, 1.0*m, 0.0*deg, 90.0*deg);
      // SBT test case e
      // testVolume = new G4Torus("test_torus1",
      //			        0.399*m, 0.4*m, 1.0*m, 0.0*deg, 90.0*deg);
      fprintf(stderr,"OK defining a Torus1 \n");
      break;
 	
    case TORUS2: 
      testVolume = new G4Torus("test_torus2",
			       0.8*m, 1.4*m, 1.8*m, startPhi, deltaPhi);
      fprintf(stderr,"OK defining a Torus2 \n");
      break;
  
    case TRAP:
      testVolume = new G4Trap ("test_trap",
			       1.0*m, 0.0, pi,
			       2.4*m,1.0*m,2.0*m, 0.0,
			       2.4*m,1.0*m,2.0*m, pi);
      fprintf(stderr,"OK defining a Trap \n");
      break;
  
    case TRD:
      testVolume = new G4Trd("test_trd",
			   0.2*m, 0.8*m, 0.8*m, 1.2*m, 4*m) ;
      fprintf(stderr,"OK defining a Trd \n");
      break;
	
    case TUBS:
      testVolume = new G4Tubs( "test_tubs", 1.0*m, 1.2*m, 1*m, startPhi, deltaPhi );

      // SBT test - case c
      // testVolume = new G4Tubs( "test_tubs", 0.0*m, 1.0*m, 1.0*m, 0.0*deg, 90.0*deg );
      // SBT test - case e
      // testVolume = new G4Tubs( "test_tubs", 0.00999*m, 0.01001*m, 1*m, 10.0*deg, 260.0*deg );
      fprintf(stderr,"OK defining a Tubs \n");
      break;
      
  //
  // specific solids
  //
  
      
    case ELLIPS:
      testVolume = new G4Ellipsoid( "test_ellipsoid",
                                     0.5*m, 0.8*m, 1.0*m, -0.4*m, 0.8*m );
      // SBT test - case a
      // testVolume = new G4Ellipsoid( "test_ellipsoid",
      //                                1.0*m, 1.0*m, 1.0*m, 0.0*m, 0.0*m );
      // SBT test - case b
      // testVolume = new G4Ellipsoid( "test_ellipsoid",
      //                                0.5*m, 0.8*m, 1.0*m, 0.0*m, 0.0*m );
      // SBT test - case c
      // testVolume = new G4Ellipsoid( "test_ellipsoid",
      //                                0.5*m, 0.8*m, 1.0*m, -0.4*m, 10.0*m );
      // SBT test - case d
      // testVolume = new G4Ellipsoid( "test_ellipsoid",
      //                                0.5*m, 0.8*m, 1.0*m, -10.0*m, 0.8*m );
      // SBT test - case e
      // testVolume = new G4Ellipsoid( "test_ellipsoid",
      //                                0.5*m, 0.8*m, 1.0*m, -0.4*m, 0.8*m );
     break;

    case ELCONE:
      testVolume = new G4EllipticalCone( "test_elcone", 0.3, 0.6, 0.75*m, 0.25*m);
      // SBT test - case a
      // testVolume = new G4EllipticalCone( "test_elcone", 0.3, 0.3, 0.75*m, 0.25*m);
      // SBT test - case b
      // testVolume = new G4EllipticalCone( "test_elcone", 0.3, 0.3, 0.75*m, 0.75*m);
      // SBT test - case c
      // testVolume = new G4EllipticalCone( "test_elcone", 0.3, 0.6, 0.75*m, 0.25*m);
      // SBT test - case d
      //testVolume = new G4EllipticalCone( "test_elcone", 0.3, 0.6, 0.75*m, 0.75*m);
     break;

    case ELTUBE:
      testVolume = new G4EllipticalTube( "test_eltube", 0.4*m, 0.8*m, 1.0*m );
      // SBT test - case a
      //testVolume = new G4EllipticalTube( "test_eltube", 1.0*m, 1.0*m, 1.0*m );
      break;
  
    case EXTRU1:
      testVolume = CreateExtrudedSolid1();
      break;
  
    case EXTRU2:
      testVolume = CreateExtrudedSolid2();
      break;
  
    case EXTRU3:
      testVolume = CreateExtrudedSolid3();
      break;
  
    case EXTRU4:
      testVolume = CreateExtrudedSolid4();
      break;
  
    case HYPE:
      testVolume = new G4Hype("test_hype", 0.2*m, 0.3*m, 0.7*rad, 0.7*rad, 0.5*m);
      // SBT test - case c
      testVolume = new G4Hype("test_hype", 0.5*m, 1.0*m, 2.0*rad, 2.0*rad, 1.0*m);
      break;
  
    case PCON: {
      G4double  z_values[2] = { -1.0*m, 1.0*m };
      G4double  rmin[2]     = {  1.0*m, 1.2*m };
      G4double  rmax[2]     = {  1.2*m, 1.4*m };
      testVolume = new G4Polycone( "test_pcon",
                                   startPhi, deltaPhi, 2, z_values, rmin, rmax );
/*
      // SBT test case d
      G4double rv[17] = { 0.0*m, 0.2*m, 0.3*m, 0.32*m, 0.32*m, 0.4*m, 0.4*m, 0.5*m, 0.5*m, 0.8*m, 
                          0.8*m, 0.9*m, 0.9*m, 0.8*m, 0.8*m, 0.3*m, 0.0*m };
      G4double zv[17] = { -0.5*m, -0.5*m, -1.1*m, -1.1*m, -0.4*m, -0.4*m, -1.0*m, -1.0*m, -0.4*m,
                          -1.0*m, 0.0*m, 0.0*m, 0.2*m, 0.2*m, 1.0*m, 0.0*m, 1.0 };
      testVolume = new G4Polycone ("test_pcon4", 0.0*deg, 90.0*deg, 17, rv, zv );

      // SBT test case e
      G4double rv[17] = { 0.0*m, 0.2*m, 0.3*m, 0.32*m, 0.32*m, 0.4*m, 0.4*m, 0.5*m, 0.5*m, 0.8*m, 
                          0.8*m, 0.9*m, 0.9*m, 0.8*m, 0.8*m, 0.3*m, 0.0*m };
      G4double zv[17] = { -0.5*m, -0.5*m, -1.1*m, -1.1*m, -0.4*m, -0.4*m, -1.0*m, -1.0*m, -0.4*m, 
                          -1.0*m, 0.0*m, 0.0*m, 0.2*m, 0.2*m, 1.0*m, 0.0*m, 1.0*m };
      testVolume = new G4Polycone ("test_pcon4", -1.0*deg, 2.0*deg, 17, rv, zv );

      // SBT test case f
      G4double rv[10] = { 0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.1*m, 0.9*m, 0.0*m, 0.0*m, 0.4*m, 0.5*m };
      G4double zv[10] = { -1.0*m, 0.0*m, 0.0*m, 0.8*m, 1.0*m, 1.0*m, 0.8*m, 0.0*m, 0.0*m, -1.0*m };
      testVolume = new G4Polycone ("test_pcon4", 10.0*deg, 250.0*deg, 10, rv, zv );
*/
      }
      break;
  
    case PCON2: {
      G4double  z_values[5] = { -1.0*m, 0.0*m, 0.0*m, 0.8*m, 1.0*m };
      G4double  rmin[5]     = {  0.5*m, 0.4*m, 0.0*m, 0.0*m, 0.9*m };
      G4double  rmax[5]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polycone( "test_pcon2",
                                   startPhi, deltaPhi, 5, z_values, rmin, rmax );
      }                                   
      break;
		
    case PCON3: {
      G4double  z_values[8] = { -1.0*m, -0.5*m, -0.5*m, -1.0*m, -1.0*m,  0.7*m,  0.7*m,  1.0*m };
      G4double  rmin[8]     = {  0.6*m,  0.6*m,  0.5*m,  0.5*m,  0.4*m,  0.4*m,  0.4*m,  0.0*m };
      G4double  rmax[8]     = {  0.7*m,  0.7*m,  0.8*m,  0.9*m,  1.0*m,  1.0*m,  0.5*m,  0.5*m };
      testVolume = new G4Polycone( "test_pcon3", 
                                   startPhi, deltaPhi, 8, z_values, rmin, rmax );
      }                                   
      break;
		
    case PCON4: {
      G4double RMINVec[8];
      RMINVec[0] = 30*cm;
      RMINVec[1] = 30*cm;
      RMINVec[2] =  0*cm;
      RMINVec[3] =  0*cm;
      RMINVec[4] =  0*cm; 
      RMINVec[5] =  0*cm;
      RMINVec[6] = 40*cm;
      RMINVec[7] = 40*cm;  

      G4double RMAXVec[8];
      RMAXVec[0] = 70*cm;
      RMAXVec[1] = 70*cm;
      RMAXVec[2] = 70*cm;
      RMAXVec[3] = 40*cm;
      RMAXVec[4] = 40*cm;
      RMAXVec[5] = 80*cm;
      RMAXVec[6] = 80*cm;
      RMAXVec[7] = 60*cm; 

      G4double Z_Values[8];
      Z_Values[0] =-20*cm;
      Z_Values[1] =-10*cm;
      Z_Values[2] =-10*cm;
      Z_Values[3] =  0*cm;
      Z_Values[4] = 10*cm;
      Z_Values[5] = 20*cm;
      Z_Values[6] = 30*cm;
      Z_Values[7] = 40*cm;

      testVolume = new G4Polycone ("test_pcon4",
                                   startPhi, deltaPhi, 8, Z_Values, RMINVec, RMAXVec );
      }                                   
      break;

    case PGON: {
      G4double z_values[2] = { -1.0*m, 1.0*m };
      G4double rmin[2]     = {  0.8*m, 1.0*m };
      G4double rmax[2]     = {  1.0*m, 1.2*m };
      testVolume = new G4Polyhedra( "test_pgon",
                                  startPhi, deltaPhi, numSide, 2, z_values, rmin, rmax );
/*
      // SBT test case c
      G4double rv[17] = { 0.0*m, 0.2*m, 0.3*m, 0.32*m, 0.32*m, 0.4*m, 0.4*m, 0.5*m, 0.5*m, 0.8*m, 
                          0.8*m, 0.9*m, 0.9*m, 0.8*m, 0.8*m, 0.3*m, 0.0*m };
      G4double zv[17] = { -0.5*m, -0.5*m, -1.1*m, -1.1*m, -0.4*m, -0.4*m, -1.0*m, -1.0*m, -0.4*m, 
                          -1.0*m, 0.0*m, 0.0*m, 0.2*m, 0.2*m, 1.0*m, 0.0*m, 1.0*m };
      testVolume = new G4Polyhedra ("test_pgon", 0.0*deg, 360.0*deg, 6, 17, rv, zv );

      // SBT test case d
      G4double rv[17] = { 0.0*m, 0.2*m, 0.3*m, 0.32*m, 0.32*m, 0.4*m, 0.4*m, 0.5*m, 0.5*m, 0.8*m, 
                          0.8*m, 0.9*m, 0.9*m, 0.8*m, 0.8*m, 0.3*m, 0.0*m };
      G4double zv[17] = { -0.5*m, -0.5*m, -1.1*m, -1.1*m, -0.4*m, -0.4*m, -1.0*m, -1.0*m, -0.4*m, 
                          -1.0*m, 0.0*m, 0.0*m, 0.2*m, 0.2*m, 1.0*m, 0.0*m, 1.0*m };
      testVolume = new G4Polyhedra ("test_pgon", 0.0*deg, 90.0*deg, 2, 17, rv, zv );

      // SBT test case f
      G4double   zv[6] = { -0.6*m, 0.0*m, -1.0*m, 0.5*m, 0.5*m, 1.0*m };
      G4double rmin[6] = { 0.5*m, 0.5*m, 0.4*m, 0.4*m, 0.8*m, 0.8*m};
      G4double rmax[6] = { 0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polyhedra ("test_pgon", 0.0*deg, 270.0*deg, 6, 6, zv, rmin, rmax );
*/
      }                                  
      break;
		
    case PGON2: {
      G4double z_values[5] = { -1.0*m, 0.0*m, 0.0*m, 0.8*m, 1.0*m };
      G4double rmin[5]     = {  0.5*m, 0.4*m, 0.0*m, 0.0*m, 0.9*m };
      G4double rmax[5]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polyhedra( "test_pgon2",
                                  startPhi, deltaPhi, numSide, 5, z_values, rmin, rmax );
      }                                  
      break;

    case PGON3: {
      G4double z_values[6] = { -0.6*m, 0.0*m,-1.0*m, 0.5*m, 0.5*m, 1.0*m };
      G4double rmin[6]     = {  0.5*m, 0.5*m, 0.4*m, 0.4*m, 0.8*m, 0.8*m };
      G4double rmax[6]     = {  0.6*m, 0.6*m, 1.0*m, 1.0*m, 1.0*m, 1.1*m };
      testVolume = new G4Polyhedra( "test_pgon3",
                                    startPhi, deltaPhi, numSide, 6, z_values, rmin, rmax );
      }                                    
      break;

    case TESSEL1: {
      // Extruded solid with triangular polygon
      // (the same as test_xtru1, but redefined via tessels

      G4TessellatedSolid* tessel = new G4TessellatedSolid("test_tessel1");
      
      tessel->AddFacet(
                new G4TriangularFacet(G4ThreeVector(-0.3*m,-0.3*m,-0.3*m),
                                      G4ThreeVector( 0.0*m, 0.3*m,-0.3*m),
                                      G4ThreeVector( 0.3*m,-0.3*m,-0.3*m),
                                      ABSOLUTE));
      tessel->AddFacet(
                new G4TriangularFacet(G4ThreeVector( 0.3*m,-0.3*m, 0.3*m),
                                      G4ThreeVector( 0.0*m, 0.3*m, 0.3*m),
                                      G4ThreeVector(-0.3*m,-0.3*m, 0.3*m),
                                      ABSOLUTE));
      tessel->AddFacet(
                new G4QuadrangularFacet(G4ThreeVector( 0.0*m, 0.3*m,-0.3*m),
                                        G4ThreeVector(-0.3*m,-0.3*m,-0.3*m),
                                        G4ThreeVector(-0.3*m,-0.3*m, 0.3*m),
                                        G4ThreeVector( 0.0*m, 0.3*m, 0.3*m),
                                        ABSOLUTE));
      tessel->AddFacet(
                new G4QuadrangularFacet(G4ThreeVector( 0.3*m,-0.3*m,-0.3*m),
                                        G4ThreeVector( 0.0*m, 0.3*m,-0.3*m),
                                        G4ThreeVector( 0.0*m, 0.3*m, 0.3*m),
                                        G4ThreeVector( 0.3*m,-0.3*m, 0.3*m),
                                        ABSOLUTE));
      tessel->AddFacet(
                new G4QuadrangularFacet(G4ThreeVector(-0.3*m,-0.3*m,-0.3*m),
                                        G4ThreeVector( 0.3*m,-0.3*m,-0.3*m),
                                        G4ThreeVector( 0.3*m,-0.3*m, 0.3*m),
                                        G4ThreeVector(-0.3*m,-0.3*m, 0.3*m),
                                        ABSOLUTE));
      tessel->SetSolidClosed(true);
      testVolume = tessel;
      G4cout << *tessel << G4endl;
/*
      G4ExtrudedSolid* xtru1 = CreateExtrudedSolid1();
      testVolume = new G4TessellatedSolid(*xtru1);
      testVolume->SetName("test_tessel1"); 
      G4cout << *((G4TessellatedSolid*)testVolume) << G4endl;
*/
      }
      break;
  
    case TESSEL2: {
/*
      G4ExtrudedSolid* xtru2 = CreateExtrudedSolid2();
      testVolume = new G4TessellatedSolid(*xtru2);
      testVolume->SetName("test_tessel2"); 
      delete xtru2;
*/
      G4ExtrudedSolid* xtru2 = CreateExtrudedSolid2();
      testVolume = xtru2;
      testVolume->SetName("test_tessel2"); 
      }                   
      break;
  
    case TESSEL3: {
/*
      G4ExtrudedSolid* xtru3 = CreateExtrudedSolid3();
      testVolume = new G4TessellatedSolid(*xtru3);
      testVolume->SetName("test_tessel3"); 
      delete xtru3;
*/
      G4ExtrudedSolid* xtru3 = CreateExtrudedSolid3();
      testVolume = xtru3;
      testVolume->SetName("test_tessel3"); 
      }                     
      break;
  
    case TESSEL4: {
/*
      G4ExtrudedSolid* xtru4 = CreateExtrudedSolid4();
      testVolume = new G4TessellatedSolid(*xtru4);
      testVolume->SetName("test_tessel4"); 
      delete xtru4;
*/
      G4ExtrudedSolid* xtru4 = CreateExtrudedSolid4();
      testVolume = xtru4;
      testVolume->SetName("test_tessel4"); 
      }                     
      break;
 
    case TET:
      testVolume = new G4Tet( "test_tet", 
                              G4ThreeVector( 0.0*m,  0.0*m,  1.0*m),
                              G4ThreeVector(-1.0*m, -1.0*m, -1.0*m),
                              G4ThreeVector(+1.0*m, -1.0*m, -1.0*m),
                              G4ThreeVector( 0.0*m,  1.0*m, -1.0*m));
/*
      // SBT test case d
      testVolume = new G4Tet( "test_tet", 
                              G4ThreeVector( 0.0*m, 0.0*m, 1.73205080756887719*m),
                              G4ThreeVector( 0.0*m, 1.63299316185545207*m, -0.577350269189625842*m),
                              G4ThreeVector(-1.41421356237309515*m, -0.816496580927726034*m, -0.577350269189625842*m),
                              G4ThreeVector( 1.41421356237309515*m, -0.816496580927726034*m, -0.577350269189625842*m));
*/
      break;

    case TWBOX:
      // SBT test case a
      testVolume = new G4TwistedBox( "test_twbox", 30.0*deg, 1.0*m,  1.0*m,  1.0*m);
      break;

    case TWTRAP1:
       testVolume = new G4TwistedTrap( "test_twtrap1", 
                                      30.0*deg, 0.8*m, 1.0*m, 1.0*m, 1.0*m );
      break;

    case TWTRAP2:
      // SBT test case a
      testVolume = new G4TwistedTrap( "test_twtrap2",
                                      30.0*deg, 1268.0*mm, 0.0*deg, 0.0*deg, 
                                      295.0*mm, 1712.2*mm, 1870.29*mm, 
                                      295.0*mm, 1712.2*mm, 1870.29*mm, 0.0*deg );
      break;

    case TWTRD:
      // SBT test case c
      testVolume = new G4TwistedTrd( "test_twtrd",
                                      0.5*m, 1.5*m, 0.25*m, 1.0*m, 1.0*m, 30.0*deg );
      break;

    case TWTUBS:
      // SBT test case a
      testVolume = new G4TwistedTubs( "test_twtubs",
                                      30.0*deg, 0.8*m, 1.0*m, -1.0*m, 1.0*m, 1, 90.0*deg );
      break;


  //
  // Boolean solids
  //

    case BOOL1:
      G4Box	*outside = new G4Box( "testboolout", 1*m, 1*m, 1*m );
      G4Tubs 	*inside = new G4Tubs( "testboolin", 0.0, 0.4*m, 1*m, 0, 360*deg );
      G4Transform3D tran = G4Translate3D( 0.4*m, 0.0, 0.0 );

      testVolume = new G4SubtractionSolid( "testbool", (G4VSolid *)outside, (G4VSolid *)inside, tran );
      break;

  }

  G4LogicalVolume	  *testLog  = new G4LogicalVolume( testVolume, Vaccuum, "test_log", 0, 0, 0 );
  new G4PVPlacement( rot, G4ThreeVector(), testLog, 
												   "test", hallLog, false, 0 );

  //
  // Put some stuff in it, if we want
  //
  if (messenger->SelectedVolume() == VOXEL) {
    G4RotationMatrix	*noRot = new G4RotationMatrix();

    G4Box 		*vxBox = new G4Box( "voxel_x", 0.3*mm, 0.6*m, 0.6*m );
    G4LogicalVolume   	*vxLog1 = new G4LogicalVolume( vxBox, Vaccuum, "x1", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( -0.6*m, 0.0*m, 0.0*m ),
												  vxLog1, "testx1", testLog, false, 0 );
    G4LogicalVolume   	*vxLog2 = new G4LogicalVolume( vxBox, Vaccuum, "x2", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( -0.2*m, 0.0*m, 0.0*m ),
												  vxLog2, "testx2", testLog, false, 0 );
    G4LogicalVolume   	*vxLog3 = new G4LogicalVolume( vxBox, Vaccuum, "x3", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( +0.2*m, 0.0*m, 0.0*m ),
												  vxLog3, "testx3", testLog, false, 0 );
    G4LogicalVolume   	*vxLog4 = new G4LogicalVolume( vxBox, Vaccuum, "x4", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( +0.6*m, 0.0*m, 0.0*m ),
												  vxLog4, "testx4", testLog, false, 0 );

    G4Box 		*vyBox = new G4Box( "voxel_y", 0.8*m, 0.3*mm, 0.6*m );
    G4LogicalVolume   	*vyLog1 = new G4LogicalVolume( vyBox, Vaccuum, "y1", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, -0.8*m, 0.0*m ),
												  vyLog1, "testy1", testLog, false, 0 );
    G4LogicalVolume   	*vyLog2 = new G4LogicalVolume( vyBox, Vaccuum, "y2", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, -0.7*m, 0.0*m ),
												  vyLog2, "testy2", testLog, false, 0 );
    G4LogicalVolume   	*vyLog3 = new G4LogicalVolume( vyBox, Vaccuum, "y3", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, +0.7*m, 0.0*m ),
												  vyLog3, "testy3", testLog, false, 0 );
    G4LogicalVolume   	*vyLog4 = new G4LogicalVolume( vyBox, Vaccuum, "y4", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, +0.8*m, 0.0*m ),
												  vyLog4, "testy4", testLog, false, 0 );

    G4Box 		*vzBox = new G4Box( "voxel_z", 0.8*m, 0.8*m, 0.3*mm );
    G4LogicalVolume   	*vzLog1 = new G4LogicalVolume( vzBox, Vaccuum, "z1", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, -0.8*m ),
												  vzLog1, "testz1", testLog, false, 0 );
    G4LogicalVolume   	*vzLog2 = new G4LogicalVolume( vzBox, Vaccuum, "z2", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, -0.7*m ),
												  vzLog2, "testz2", testLog, false, 0 );
    G4LogicalVolume   	*vzLog3 = new G4LogicalVolume( vzBox, Vaccuum, "z3", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, +0.7*m ),
												  vzLog3, "testz3", testLog, false, 0 );
    G4LogicalVolume   	*vzLog4 = new G4LogicalVolume( vzBox, Vaccuum, "z4", 0, 0, 0 );
    new G4PVPlacement( noRot, G4ThreeVector( 0.0*m, 0.0*m, +0.8*m ),
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
