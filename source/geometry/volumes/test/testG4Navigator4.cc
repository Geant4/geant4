// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Navigator4.cc,v 1.3 1999-12-15 14:50:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//   Locate & Step within simple replicated geometry, both
//   with and without voxels.

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("cube",50,50,50);
    G4Box *myBox=new G4Box("large cube",25,25,25);
    G4Box *mySlab=new G4Box("slab",6.25,25,25);
    G4Box *myBox2=new G4Box("small cube",10,10,10);
    G4Box *mySlab2=new G4Box("small slab",10,2.5,10);
    G4Box *mySlab3=new G4Box("oblong",2.5,2.5,10);

    G4Box *myTargetBox=new G4Box ("small cube",2,2,2);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
    // Note: no mother pointer set
    
// Box containing x replicas centred on origin
    G4LogicalVolume *boxLog=new G4LogicalVolume(myBox,0,
						"Container",0,0,0);
    
    G4PVPlacement *boxPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					     "Container",boxLog,
					     worldPhys,false,0);
    G4LogicalVolume *repLog=new G4LogicalVolume(mySlab,0,
						"Slab",0,0,0);
    G4PVReplica *slabPhys=new G4PVReplica("Slab",
					  repLog,boxPhys,
					  kXAxis,4,12.5);

// Box containing y then x replicas centred on 0,0,35 (touches x replica
// box (above)

    G4LogicalVolume *box2Log=new G4LogicalVolume(myBox2,0,
						 "Container2",0,0,0);
    
    G4PVPlacement *box2Phys=new G4PVPlacement(0,G4ThreeVector(0,0,35),
					      "Container2",box2Log,
					      worldPhys,false,0);
    G4LogicalVolume *rep2Log=new G4LogicalVolume(mySlab2,0,
						 "Slab2",0,0,0);
    G4PVReplica *slab2Phys=new G4PVReplica("Slab2",
					   rep2Log,box2Phys,
					   kYAxis,4,5);
    G4LogicalVolume *rep3Log=new G4LogicalVolume(mySlab3,0,
						 "Slab3",0,0,0);
    G4PVReplica *slab3Phys=new G4PVReplica("Slab3",
					   rep3Log,slab2Phys,
					   kXAxis,4,5);


// Box containing y then x replicas centred (35,35,35)
// (touches x replica box at one vertex)
// Cube positioned at centre of x replicas

    G4LogicalVolume *box3Log=new G4LogicalVolume(myBox2,0,
						 "Container3",0,0,0);
    
    G4RotationMatrix *rot=new G4RotationMatrix();
    rot->rotateZ(M_PI/2);
    G4PVPlacement *box3Phys=new G4PVPlacement(rot,G4ThreeVector(35,35,35),
					      "Container3",box3Log,
					      worldPhys,false,0);
    G4LogicalVolume *rep4Log=new G4LogicalVolume(mySlab2,0,
						"Slab4",0,0,0);
    G4PVReplica *slab4Phys=new G4PVReplica("Slab4",
					   rep4Log,box3Phys,
					   kYAxis,4,5);
    G4LogicalVolume *rep5Log=new G4LogicalVolume(mySlab3,0,
						 "Slab5",0,0,0);
    G4PVReplica *slab5Phys=new G4PVReplica("Slab5",
					   rep5Log,slab4Phys,
					   kXAxis,4,5);
    G4LogicalVolume  *targetLog=new G4LogicalVolume(myTargetBox,0,
						    "Target",0,0,0);
    G4PVPlacement *targetPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
						"Target",targetLog,
						slab5Phys,false,0);

    return worldPhys;
}


//
// Test LocateGlobalPointAndSetup
//
G4bool testG4NavigatorLocate(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,0,0),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(1,0,0),0,false);
    assert(located->GetName()=="Slab");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(2,0,0));
    assert(located->GetName()=="Slab");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,0,0));
    assert(located->GetName()=="World");


    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-19.75,10,-10));
    assert(located->GetName()=="Slab");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(-1,10,-10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7.25,10,-10));
    assert(located->GetName()=="Slab");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(-1,10,-10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(5.25,10,-10));
    assert(located->GetName()=="Slab");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(-1,10,-10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(17.75,10,-10));
    assert(located->GetName()=="Slab");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(-1,10,-10)));

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7.5,7.5,45),0,false);
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(0,0,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-1.5,7.5,45),0,false);
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(1,0,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(4.5,7.5,45),0,false);
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(2,0,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(10,7.5,45),0,false);
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(2.5,0,10)));


    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7.5,-7.5,45));
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(0,0,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-1.5,-1.5,45));
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(1,1,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(4.5,4.5,45));
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(2,2,10)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(10,10,45));
    assert(located->GetName()=="Slab3");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(2.5,2.5,10)));

    // Check found to be inside daughter volumes of replicas
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(36,36,36),0,false);
    //    G4cout << "At " << G4ThreeVector(36,36,36) << G4endl << myNav;
    assert(located->GetName()=="Target");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(1.5,-1.5,1)));
    return true;
}

//
// Test ComputeStep
//
G4bool testG4NavigatorSteps(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector pos,dir,origin,xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    myNav.SetWorldVolume(pTopNode);

    pos=G4ThreeVector(1,0,0);
    dir=yHat;
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="Slab");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,25));
    assert(ApproxEqual(safety,1));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,25));
    assert(safety==0);

    pos=G4ThreeVector(1,1,-50);
    dir=zHat;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,25));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,50));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab3");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));

    Step=10;
    pos+=Step*dir;
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab3");

    dir=myHat;
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,1));
    assert(ApproxEqual(safety,1));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab3");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab3");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="World");


    pos=G4ThreeVector(36,36,36);
    dir=yHat;
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Target");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,3.5));
    assert(ApproxEqual(safety,0.5));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab5");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,0.5));
    assert(ApproxEqual(safety,0));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab5");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,0.5));
    assert(ApproxEqual(safety,0));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Target");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    //    G4cout << "At " << pos << G4endl << myNav;
    assert(located->GetName()=="Slab5");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(-2,-1.5,1)));
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,0.5));
    assert(ApproxEqual(safety,0));

    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    
    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4NavigatorLocate(myTopNode);
    testG4NavigatorSteps(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4NavigatorLocate(myTopNode);
    testG4NavigatorSteps(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore *ps=G4PhysicalVolumeStore::GetInstance();
    for (G4int i=ps->entries()-1;i>=0;i--)
      {
	// Delete any rotation matrices
	delete ps->at(i)->GetRotation();
      }
    return 0;
}




