// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Navigator5.cc,v 1.1 1999-01-08 16:31:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Locate & Step within replicated geometry (without full voxels)
// consisting of tubes (replicated in r,phi...)     P.Kent August 96

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

// Resolve template problem on some architectures.
#include "locg4templates.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4VPhysicalVolume* BuildGeometry()
{

    G4Box *worldBox= new G4Box ("cube",50,50,50);
    G4Tubs *fullTube= new G4Tubs("tube",0,10,10,0,2*M_PI);
    G4Tubs *fullTubeSlice= new G4Tubs("tube slice",0,2,10,0,2*M_PI);

    G4Tubs *hollowTube= new G4Tubs("hollow tube",2,8,10,0,2*M_PI);
    G4Tubs *hollowTubeSlice= new G4Tubs("hollow tube slice",2,4,10,0,2*M_PI);

    G4Tubs *hphiTube= new G4Tubs("phi tube",2,8,10,-M_PI*0.5,M_PI);
    G4Tubs *hphiTubeSlice= new G4Tubs("phi tube slice",2,8,10,0,M_PI*0.25);


    G4Tubs *allTube= new G4Tubs("all tube",2,10,10,0,2*M_PI);
    G4Tubs *allTubeZSlice= new G4Tubs("all tube z slice",2,10,2,0,2*M_PI);
    G4Tubs *allTubeZRSlice=new G4Tubs("all tube zr slice",2,6,2,0,2*M_PI);
    G4Tubs *allTubeZRPSlice=new G4Tubs("all tube zrp slice",2,6,
				       2,-M_PI*0.25,M_PI*0.5);


    G4LogicalVolume *worldLog=new G4LogicalVolume(worldBox,0,
						  "World",0,0,0);
    // Logical with no field,sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
    // Note: no mother pointer set
    

    // Full (solid) Tube containing r replicas
    G4LogicalVolume *fullTubeLog=new G4LogicalVolume(fullTube,
						     0,
						     "Container",
						     0,0,0);
    G4PVPlacement *fullTubePhys=new G4PVPlacement(0,G4ThreeVector(0,0,35),
						  "Container",
						  fullTubeLog,
						  worldPhys,false,0);
    G4LogicalVolume *fullTubeSliceLog=new G4LogicalVolume(fullTubeSlice,
							  0,
							  "Slice",
							  0,0,0);
    G4PVReplica *fullTubeSlicePhys=new G4PVReplica("TubeSlice",
					  fullTubeSliceLog,fullTubePhys,
					  kRho,5,2,0);


    // hollow Tube containing r replicas
    G4LogicalVolume *hollowTubeLog=new G4LogicalVolume(hollowTube,
						       0,
						       "Container",
						       0,0,0);
    G4PVPlacement *hollowTubePhys=new G4PVPlacement(0,G4ThreeVector(0,0,-35),
						    "Container",
						    hollowTubeLog,
						    worldPhys,false,0);

    G4LogicalVolume *hollowTubeSliceLog=new G4LogicalVolume(hollowTubeSlice,
							    0,
							    "Slice",
							    0,0,0);
    G4PVReplica *hollowTubeSlicePhys=new G4PVReplica("HollowTubeSlice",
					  hollowTubeSliceLog,hollowTubePhys,
					  kRho,3,2,2);


    // Hollow Tube containing phi replicas
    G4LogicalVolume *hphiTubeLog=new G4LogicalVolume(hphiTube,
						     0,
						     "Container",
						     0,0,0);
    G4PVPlacement *hphiTubePhys=new G4PVPlacement(0,G4ThreeVector(35,0,0),
						  "Container",
						  hphiTubeLog,
						  worldPhys,false,0);

    G4LogicalVolume *hphiTubeSliceLog=new G4LogicalVolume(hphiTubeSlice,
							  0,
							  "Slice",
							  0,0,0);
    G4PVReplica *hphiTubeSlicePhys=new G4PVReplica("hphiTubeSlice",
					  hphiTubeSliceLog,hphiTubePhys,
					  kPhi,4,M_PI*0.25,-M_PI*0.5);

    // Tube replicated along z, then r then phi (!)

    G4LogicalVolume *allTubeLog=new G4LogicalVolume(allTube,
						    0,
						    "Container",
						    0,0,0);
    G4PVPlacement *allTubePhys=new G4PVPlacement(0,G4ThreeVector(0,35,0),
						 "Container",
						 allTubeLog,
						 worldPhys,false,0);

    G4LogicalVolume *allTubeZSliceLog=new G4LogicalVolume(allTubeZSlice,
							  0,
							  "Slice",
							  0,0,0);
    G4PVReplica *allTubeZSlicePhys=new G4PVReplica("allTubeZSlice",
						   allTubeZSliceLog,
						   allTubePhys,
						   kZAxis,5,4);

    G4LogicalVolume *allTubeZRSliceLog=new G4LogicalVolume(allTubeZRSlice,
							   0,
							   "Slice",
							   0,0,0);
    G4PVReplica *allTubeZRSlicePhys=new G4PVReplica("allTubeZRSlice",
						    allTubeZRSliceLog,
						    allTubeZSlicePhys,
						    kRho,2,4,2);

    G4LogicalVolume *allTubeZRPSliceLog=new G4LogicalVolume(allTubeZRPSlice,
							    0,
							    "Slice",
							    0,0,0);
    G4PVReplica *allTubeZRPSlicePhys=new G4PVReplica("allTubeZRPSlice",
						     allTubeZRPSliceLog,
						     allTubeZRSlicePhys,
						     kPhi,4,M_PI*0.5,-M_PI*0.25);


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
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(50,50,50),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,35),0,false);
    assert(located->GetName()=="TubeSlice");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(0,0,0)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-35),0,false);
    assert(located->GetName()=="World");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(4,0,-35),0,false);
    assert(located->GetName()=="HollowTubeSlice");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(4,0,0)));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(37,0.5,0),0,false);
    assert(located->GetName()=="hphiTubeSlice");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(3,35,0),0,false);
    assert(located->GetName()=="allTubeZRPSlice");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),
		       G4ThreeVector(3,0,0)));
    
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

    pos=G4ThreeVector(0,0,49);
    dir=mzHat;
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,1));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(safety==0);
    Step=10;
    pos+=Step*dir;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,10));
    assert(ApproxEqual(safety,2));

    pos=G4ThreeVector(-10,-10,35);
    dir=(xHat+yHat).unit();
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,sqrt(200.)-10.));
    assert(ApproxEqual(safety,sqrt(200.)-10.));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="TubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
  


    pos=G4ThreeVector(10,10,-35);
    dir=(mxHat+myHat).unit();
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,sqrt(200.)-8.));
    assert(ApproxEqual(safety,sqrt(200.)-8.));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="HollowTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
 


    pos=G4ThreeVector(39,-10,0);
    dir=yHat;
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,10.-sqrt(48.)));
    assert(ApproxEqual(safety,sqrt(116.)-8.));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="hphiTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,sqrt(48.)-4.));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="hphiTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="hphiTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="hphiTubeSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,sqrt(48.)-4.));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");


    pos=G4ThreeVector(12,35,0);
    dir=mxHat;
    physStep=kInfinity;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,2));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="allTubeZRPSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="allTubeZRPSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    Step=2;
    pos+=Step*dir;
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="allTubeZRPSlice");
    dir=zHat;
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,2));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="allTubeZRPSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="allTubeZRPSlice");
    Step=myNav.ComputeStep(pos,dir,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    pos+=Step*dir;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(pos);
    assert(located->GetName()=="World");
  
    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4NavigatorLocate(myTopNode);
    testG4NavigatorSteps(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}




