// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Navigator3.cc,v 1.2 1999-11-24 21:13:53 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//   Locate & Step within simple rotated boxlike geometry, both
//   with and without voxels.

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("cube",50,50,50);
    G4Box *myCuboid=new G4Box("cuboid",5,10,15);
    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
    // Note: no mother pointer set
    
    G4LogicalVolume *cubLog=new G4LogicalVolume(myCuboid,0,
						"Crystal Box",0,0,0);
    
    G4RotationMatrix *rot1=new G4RotationMatrix();
    rot1->rotateZ(M_PI*0.5);
    G4PVPlacement *cubPhys1=new G4PVPlacement(rot1,G4ThreeVector(0,0,10),
					      "Target 1",cubLog,
					      worldPhys,false,0);

    G4RotationMatrix *rot2=new G4RotationMatrix();
    rot2->rotateX(M_PI*0.5);
    G4PVPlacement *cubPhys2=new G4PVPlacement(rot2,G4ThreeVector(-30,0,10),
					      "Target 2",cubLog,
					      worldPhys,false,0);
    G4RotationMatrix *rot3=new G4RotationMatrix();
    rot3->rotateY(M_PI*0.5);
    G4PVPlacement *cubPhys3=new G4PVPlacement(rot3,G4ThreeVector(30,0,10),
					      "Target 3",cubLog,
					      worldPhys,false,0);
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
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,10),0,false);
    assert(located->GetName()=="Target 1");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(10,0,10),0,false);
    assert(located->GetName()=="Target 1");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,10,10),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,5,10),0,false);
    assert(located->GetName()=="Target 1");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-30,0,10),0,false);
    assert(located->GetName()=="Target 2");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-30,10,25),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-35,15,20),0,false);
    assert(located->GetName()=="Target 2");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25,-15,0),0,false);
    assert(located->GetName()=="Target 2");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,0,10),0,false);
    assert(located->GetName()=="Target 3");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(30,15,15),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,10,15),0,false);
    assert(located->GetName()=="Target 3");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(45,-10,5),0,false);
    assert(located->GetName()=="Target 3");

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

    pos=G4ThreeVector(0,0,10);
    dir=xHat;
    located=myNav.LocateGlobalPointAndSetup(pos,0,false);
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,10));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();

    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,30));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,5));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located==0);
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




