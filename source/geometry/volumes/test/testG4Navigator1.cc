// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Navigator1.cc,v 1.2 1999-11-24 21:13:51 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

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

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Sample Paramterisation
class G4LinScale : public G4VPVParameterisation
{
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,(n-1)*15,0));
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int n,
				 const G4VPhysicalVolume* pRep) const
  {
    pBox.SetXHalfLength(10);
    pBox.SetYHalfLength(5+n);
    pBox.SetZHalfLength(5+n);
  }
};
G4LinScale myParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("cuboid",25,25,20);
    G4Box *myBox=new G4Box("cube",10,10,10);
    G4Box *mySlab= new G4Box("slab",10,25,10);

    G4Box *myVariableBox=new G4Box("Variable Box",10,5,5);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume *boxLog=new G4LogicalVolume(myBox,0,
						"Crystal Box",0,0,0);
    G4LogicalVolume *slabLog=new G4LogicalVolume(mySlab,0,
						 "Crystal Slab",0,0,0);
    G4PVPlacement *offMXYPhys=new G4PVPlacement(0,G4ThreeVector(-15,15,-10),
						"Target 1",boxLog,
						worldPhys,false,0);
    G4PVPlacement *offMXMYPhys=new G4PVPlacement(0,G4ThreeVector(-15,-15,-10),
						 "Target 2",boxLog,
						 worldPhys,false,0);

    G4PVPlacement *offYPhys=new G4PVPlacement(0,G4ThreeVector(15,0,-10),
					      "Target 3",slabLog,
					      worldPhys,false,0);

    G4PVPlacement *offYZPhys=new G4PVPlacement(0,G4ThreeVector(0,15,10),
					       "Target 4",boxLog,
					       worldPhys,false,0);
    G4PVPlacement *offMYZPhys=new G4PVPlacement(0,G4ThreeVector(0,-15,10),
						"Target 5",boxLog,
						worldPhys,false,0);


    G4LogicalVolume *variLog=new G4LogicalVolume(myVariableBox,0,
					         "Variable Blocks",0,0,0);
    G4PVParameterised *paramPhys=new G4PVParameterised("Vari' Blocks",
						       variLog,
						       offYPhys,
						       kYAxis,
						       3,
						       &myParam);
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

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check relative search that causes backup one level and then search down:
// Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10),0,false);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-15,20));
    assert(located->GetName()=="Target 5");
    assert(ApproxEqual(myNav.GetCurrentLocalCoordinate(),G4ThreeVector(0,0,10)));
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check parameterised volumes

// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-16));
    assert(located->GetName()=="Target 3");

// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-17));
    assert(located->GetName()=="Target 3");

// Replication 2
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-18));
    assert(located->GetName()=="Target 3");

    return true;
}


//
// Test Stepping
//
G4bool testG4Navigator2(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    myNav.SetWorldVolume(pTopNode);
  
//
// Test location & Step computation
//  
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mxHat,physStep,safety);
    assert(ApproxEqual(Step,25));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),xHat,physStep,safety);
    assert(ApproxEqual(Step,5));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(5,0,-10),0,true);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),zHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mzHat,physStep,safety);
    assert(ApproxEqual(Step,10));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);


//
// Test stepping through common boundaries
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,-20));
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,-20),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,0));
    assert(located->GetName()=="Target 4");
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,0),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,20));
    assert(!located);

//
// Test mother limited Step
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25,0,10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-25,0,10),xHat,physStep,safety);
    assert(ApproxEqual(Step,50));
    assert(ApproxEqual(safety,0));

//
// Test stepping through parameterised volumes
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-25,-10),0,false);
    assert(located->GetName()=="Target 3");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(15,-25,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-20,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-20,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,10));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-10,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,-10,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-6,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,12));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,6,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,8,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,8,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,14));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,22,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,22,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,3));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,25,-10));
    assert(!located);

    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}




