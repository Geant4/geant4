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
// $Id: testExitNormalNav.cc,v 1.2 2002-10-22 12:40:30 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   Testing the product of Exit Normal of the Navigator for
//     simple hierarchial geometry.  
//      ( replicas, parameterised volumes currently not included )
//  
// First version:  J. Apostolakis,  18th June 2002

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
// #include "G4PVParameterised.hh"
// #include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Build simple geometry:
//   6 small cubes inside a slab (all G4Boxes) 
//   3 slabs are positioned inside the world cuboid

G4VPhysicalVolume* BuildGeometry()
{

    // Rotations in X
  G4RotationMatrix *prot90d_X, *prot180d_X, *prot270d_X;
    prot90d_X =  new G4RotationMatrix();
    prot180d_X = new G4RotationMatrix();
    prot270d_X = new G4RotationMatrix();
    prot90d_X->rotateX(M_PI*0.5);
    prot180d_X->rotateX(M_PI);
    prot270d_X->rotateX(M_PI*1.5);

    // Rotations in Y
    G4RotationMatrix *prot90d_Y, *prot180d_Y, *prot270d_Y;
    prot90d_Y =  new G4RotationMatrix();
    prot180d_Y = new G4RotationMatrix();
    prot270d_Y = new G4RotationMatrix();
    prot90d_Y->rotateY(M_PI*0.5);
    prot180d_Y->rotateY(M_PI);
    prot270d_Y->rotateY(M_PI*1.5);

    // Rotations in Z
    G4RotationMatrix *prot90d_Z, *prot180d_Z, *prot270d_Z;
    prot90d_Z =  new G4RotationMatrix();
    prot180d_Z = new G4RotationMatrix();
    prot270d_Z = new G4RotationMatrix();
    prot90d_Z->rotateZ(M_PI*0.5);
    prot180d_Z->rotateZ(M_PI);
    prot270d_Z->rotateZ(-M_PI*0.5);

    // Solids
    G4Box *myBigBox=  new G4Box("BigBox-World",200.*cm,200.*cm,200.*cm);
    G4Box *Slab=      new G4Box("slab",17.5*cm,10.*cm,7.5*cm);
    G4Box *inCube10= new G4Box("Cube ten",5.*cm,5.*cm,5.*cm);
    G4Box *smallCube= new G4Box("Small cube", 0.5*cm, 0.5*cm, 0.5*cm);

    // World
    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "WorldLV",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "WorldPV",worldLog,
					       0,false,0);
				// Note: no mother pointer set

    // Slab volume
    G4LogicalVolume *slabLog=new G4LogicalVolume(Slab, 0, "Slab-logV");    

    // Inner volume
    G4LogicalVolume *boxLog=
            new G4LogicalVolume(inCube10, 0, "Cube10-logV");

    // Smallest cube volume
    G4LogicalVolume *smallLog=
            new G4LogicalVolume(smallCube, 0, "smallCube1-lV");

    // Place small cubes inside Cube ten
    new G4PVPlacement(prot90d_Y,
		      G4ThreeVector( -4.5*cm, 0.0, 0.0),
		      smallLog,
		      "smallBackX",
		      boxLog,
		      false,
		      0);    
    new G4PVPlacement(prot180d_Y,
		      G4ThreeVector(  4.5*cm, 0.0, 0.0),
		      smallLog,
		      "smallFrontX",
		      boxLog,
		      false,
		      1);    
    new G4PVPlacement(prot90d_Z,
		      G4ThreeVector( 0.0, -4.5*cm, 0.0),
		      smallLog,
		      "smallBackY",
		      boxLog,
		      false,
		      2);    
    new G4PVPlacement(prot90d_X,
		      G4ThreeVector( 0.0,  4.5*cm, 0.0),
		      smallLog,
		      "smallFrontY",
		      boxLog,
		      false,
		      3);    

    // Fill the slab with inner volumes
    //
    
    G4ThreeVector    centerPositionFirst(12.5*cm,-5*cm,0.0);
    new G4PVPlacement(prot90d_Y,
		      centerPositionFirst,
		      boxLog,
		      "Lower Front",
		      slabLog,
		      false,
		      0);

    G4ThreeVector    centerPositionSecond(12.5*cm, 5*cm,0.0);
    new G4PVPlacement(prot180d_Z,
		      centerPositionSecond,
		      boxLog,
		      "Upper Front",
		      slabLog,
		      false,
		      1);

    G4ThreeVector    centerPositionThird(-12.5*cm, 5*cm,0.0);
    new G4PVPlacement(prot270d_X,
		      centerPositionThird,
		      boxLog,
		      "Upper Back",
		      slabLog,
		      false,
		      2);


    G4ThreeVector    centerPositionFourth(-12.0*cm, -5*cm,0.0);
    new G4PVPlacement(prot180d_Z,
		      centerPositionFourth,
		      boxLog,
		      "Lower Back",
		      slabLog,
		      false,
		      3);


    G4ThreeVector    centerPositionFifth(-2.5*cm, 5*cm,0.0);
    new G4PVPlacement(prot90d_Y,
		      centerPositionFifth,
		      boxLog,
		      "Upper Mid-Back",
		      slabLog,
		      false,
		      4);


    G4ThreeVector    centerPositionSixth( 2.5*cm, -5*cm,0.0);
    new G4PVPlacement(prot90d_X,
		      centerPositionSixth,
		      boxLog,
		      "Lower Mid-Front",
		      slabLog,
		      false,
		      5);

    // Placement of Slabs in World Volume
    // 
    G4ThreeVector    slabPositionOne( -27.5*cm, 0.0,0.0);
    new G4PVPlacement(prot180d_Y,
		      slabPositionOne,
		      "Back-Slab1",
		      slabLog,
		      worldPhys,
		      false,
		      1);

    G4ThreeVector    slabPositionTwo( 0.0, 0.0, 0.0);
    new G4PVPlacement(prot90d_Z,
		      slabPositionTwo,
		      "Upright-Middle-Slab2",
		      slabLog,
		      worldPhys,
		      false,
		      2);


    G4ThreeVector    slabPositionThree( 27.5*cm, 0.0, 0.0);
    new G4PVPlacement(prot180d_Z,
		      slabPositionThree,
		      "Front-Slab3",
		      slabLog,
		      worldPhys,
		      false,
		      3);

    return worldPhys;
}

//
// Test LocateGlobalPointAndSetup
//
G4bool testG4Navigator1(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(1000*cm,0,0),0,false));

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(1800.0*mm,0,0),0,false);
    assert(located->GetName()=="WorldPV");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(1000.*cm,0,0)));

    return true;

    // Can add more location checks here, like the old ones below.

// Check relative search that causes backup one level and then search down:
// Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15.*cm,0,-5.*cm),0,false);
    assert(located->GetName()=="Upper Front");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-15.*cm,20.*cm));
    assert(located->GetName()=="Target 5");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),G4ThreeVector(0,0,10)));
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

    return true;
}

int verbose= 1;
//
// Test Stepping
//
G4bool testExitNormal(G4VPhysicalVolume *pTopNode,
                      G4ThreeVector     initialPoint, 
		      G4ThreeVector     direction,
		      G4ThreeVector     expectedExitNorm)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    
    myNav.SetWorldVolume(pTopNode);
//
// Test location & Step computation
//  
    G4ThreeVector initPoint(initialPoint), newPoint(0,0,0);
    // G4ThreeVector direction= xHat;
    G4bool valid;

    if( verbose ){
      G4cout << "Initial step " << G4endl;
      G4cout << "-Initial Point = "  << initPoint    << G4endl;
    }

    located=myNav.LocateGlobalPointAndSetup(initPoint);
    assert(located->GetName()=="WorldPV");
    if( verbose )
      G4cout << "-Located: Location before is " << located->GetName() << G4endl; 

    physStep=kInfinity;
    Step=myNav.ComputeStep(initPoint, direction, physStep, safety);
    if( verbose ){
      G4cout << "-Moved: Step was = " << Step << " expected " << 5.0 * cm << G4endl;
      G4cout << "  safety= " << safety << G4endl;
    }
    assert(ApproxEqual(Step,5.0*cm));
    // assert(ApproxEqual(safety,50.0));
    assert(safety>=0.0);
    assert(safety<=50.0);

    newPoint= initPoint + Step * direction; 
    G4ThreeVector localNormal = myNav.GetLocalExitNormal(&valid); 
    assert(valid);
    G4ThreeVector globalNormal = myNav.GetLocalToGlobalTransform().TransformAxis(localNormal);
    assert( globalNormal == expectedExitNorm );

    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(initPoint);
    assert(located->GetName()!="WorldPV");
    if( verbose )
      G4cout << "-Located: Location is " << located->GetName() << G4endl; 

    // Next Steps
    G4int istep;
    for ( istep=0; istep < 15; istep++ ){

       initPoint= newPoint;

       if( verbose ){
	 G4cout << "Sub step " << istep << G4endl;
	 G4cout << "-Initial Point = "  << initPoint    << G4endl;
	 G4cout << "-Location before is " << located->GetName() << G4endl; 
       }

       physStep=kInfinity;
       Step=myNav.ComputeStep(initPoint, direction, physStep, safety);
       if( verbose )
	 G4cout << "-Moved: Step was = " << Step << G4endl;
       assert( Step <= 10.0*cm);
       assert(ApproxEqual(safety,0.0));
       assert(safety>=0);

       newPoint= initPoint + Step * direction; 
       G4ThreeVector localNormal = myNav.GetLocalExitNormal(&valid); 
       assert(valid);
       G4ThreeVector globalNormal = myNav.GetLocalToGlobalTransform().TransformAxis(localNormal);
      
       if( 0 ) { // globalNormal != G4ThreeVector(1.0,0.0,0.0) ){
         G4cout << " **Problem** with pre-relocation normals: " << G4endl;
	 G4cout << "   *Point      = "  << newPoint    << G4endl;
	 G4cout << "   *localNorm  = "  << localNormal << G4endl;
	 G4cout << "   *globalNorm = " << globalNormal << G4endl;
       }
       //  assert( globalNormal == G4ThreeVector(1.0,0.0,0.0) );

       myNav.SetGeometricallyLimitedStep();
       located=myNav.LocateGlobalPointAndSetup(newPoint);
       // assert(located->GetName()!="WorldPV");
       if( verbose )
	 G4cout << "-Located: Location after is " << located->GetName() << G4endl; 

       localNormal = myNav.GetLocalExitNormal(&valid); 
       assert(valid);
       globalNormal = myNav.GetLocalToGlobalTransform().TransformAxis(localNormal);
       if( verbose ) {
	 G4cout << "Post-relocation normals: " << G4endl;
	 G4cout << " Point      = "  << newPoint    << G4endl;
	 G4cout << " Location after is " << located->GetName() << G4endl; 
	 G4cout << " localNorm  = "  << localNormal << G4endl;
	 G4cout << " globalNorm = " << globalNormal << G4endl;
       }
       assert( ApproxEqual( globalNormal, expectedExitNorm ) );

    }

    return true; 
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    const G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    const G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);

    G4ThreeVector initPointMinusX(-50.0*cm,0.01*cm,0.);
    G4ThreeVector initPointPluxX(50.0*cm, -0.01*cm,0.);

    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4Navigator1(myTopNode);
    testExitNormal(myTopNode, initPointMinusX, xHat,  xHat);
    testExitNormal(myTopNode, initPointPluxX, mxHat, mxHat);
    testExitNormal(myTopNode, G4ThreeVector(-50.0*cm,2.0*cm,0.0),  xHat, xHat);
    testExitNormal(myTopNode, G4ThreeVector(-50.0*cm,-2.0*cm,0.0), xHat, xHat);

// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4Navigator1(myTopNode);
    testExitNormal(myTopNode, initPointMinusX, xHat, xHat);
    testExitNormal(myTopNode, initPointPluxX, mxHat, mxHat);
    testExitNormal(myTopNode, G4ThreeVector(-50.0*cm,2.0*cm,0.0),  xHat, xHat);
    testExitNormal(myTopNode, G4ThreeVector(-50.0*cm,-2.0*cm,0.0), xHat, xHat);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}




