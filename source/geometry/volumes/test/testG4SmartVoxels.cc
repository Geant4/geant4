// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4SmartVoxels.cc,v 1.2 1999-12-15 14:50:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// testG4SmartVoxels
//             Ensure asserts are compiled in


#include <assert.h>

#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4Box.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"

#include "G4SmartVoxelNode.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4bool testG4SmartVoxelNodes()
{
    G4SmartVoxelNode tNode(1);
    assert(tNode.GetNoContained()==0);
    assert(tNode.GetMinEquivalentSliceNo()==1);
    assert(tNode.GetMaxEquivalentSliceNo()==1);
    tNode.SetMinEquivalentSliceNo(0);
    tNode.SetMaxEquivalentSliceNo(2);
    assert(tNode.GetMinEquivalentSliceNo()==0);
    assert(tNode.GetMaxEquivalentSliceNo()==2);
    return true;

}

// Test of geometry close and hence voxel Build
// 2 small cubes are positioned inside a larger cuboid
//
// [Navigation links [logi/phys volumes] checked by testG4Volumes]
G4bool testG4SmartVoxels()
{

    G4Box myBigBox("cuboid",25,25,20);
    G4Box myBox("cube",10,10,10);
    G4Box mySlab("slab",10,25,10);

    G4LogicalVolume myDetectorLog(&myBigBox,0,
				  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement myDetectorPhys(0,G4ThreeVector(0,0,0),
				 "World",&myDetectorLog,
				 0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume myDaughter1Log(&myBox,0,"Crystal Box",0,0,0);
    G4LogicalVolume myDaughter2Log(&mySlab,0,"Crystal Slab",0,0,0);
    G4PVPlacement offMXYPhys(0,G4ThreeVector(-15,15,-10),
			     "Target 1",&myDaughter1Log,
			     &myDetectorPhys,false,0);
    G4PVPlacement offMXMYPhys(0,G4ThreeVector(-15,-15,-10),
			      "Target 2",&myDaughter1Log,
			      &myDetectorPhys,false,0);

    G4PVPlacement offYPhys(0,G4ThreeVector(15,0,-10),
			   "Target 3",&myDaughter2Log,
			   &myDetectorPhys,false,0);

    G4PVPlacement offYZPhys(0,G4ThreeVector(0,15,10),
			    "Target 4",&myDaughter1Log,
			    &myDetectorPhys,false,0);
    G4PVPlacement offMYZPhys(0,G4ThreeVector(0,-15,10),
			     "Target 5",&myDaughter1Log,
			     &myDetectorPhys,false,0);
    
// Close geometry and check voxels constructed
    assert(G4GeometryManager::GetInstance()->CloseGeometry());

// Check we have some voxels
    assert(myDetectorLog.GetVoxelHeader() != 0);
    if (kMinVoxelVolumesLevel1==0)
	{
	    assert(myDaughter1Log.GetVoxelHeader() != 0);
	    assert(myDaughter2Log.GetVoxelHeader() != 0);
	}
    else
	{
	    assert(!myDaughter1Log.GetVoxelHeader());
	    assert(!myDaughter2Log.GetVoxelHeader());
	}
//
//
//


// Open geometry freeing voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
// Check voxels destroyed
    assert(!myDetectorLog.GetVoxelHeader());
    assert(!myDaughter1Log.GetVoxelHeader());
    assert(!myDaughter2Log.GetVoxelHeader());

    return true;
}

// Test voxel Build for replicas
G4bool testG4ReplicaVoxels()
{

    G4Box myBigBox("cuboid",100,100,100);
    G4Box myBox("cube",100,100,50);
    G4Box mySlab("slab",100,100,10);

    G4LogicalVolume myDetectorLog(&myBigBox,0,
				  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement myDetectorPhys(0,G4ThreeVector(0,0,0),
				 "World",&myDetectorLog,
				 0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume myDaughter1Log(&myBox,0,"Crate",0,0,0);
    G4LogicalVolume myDaughter2Log(&mySlab,0,"Slab",0,0,0);
    G4PVPlacement myPhys1(0,G4ThreeVector(0,0,0),
			  "Target 1",&myDaughter1Log,
			  &myDetectorPhys,false,0);

    G4PVReplica myPhysRep1("Replicated Slabs",
                           &myDaughter2Log,
                           &myPhys1,
			   kZAxis,5,20);


// Close geometry and check voxels constructed
    assert(G4GeometryManager::GetInstance()->CloseGeometry());

// Check replicated voxels
    assert(myDaughter1Log.GetVoxelHeader() != 0);
    assert(!myDaughter2Log.GetVoxelHeader());


    G4SmartVoxelHeader *vHead=myDaughter1Log.GetVoxelHeader();
    G4SmartVoxelProxy *vProxy;
    G4SmartVoxelNode *vNode;
    assert(vHead->GetAxis()==kZAxis);
    assert(vHead->GetMinExtent()==-50);
    assert(vHead->GetMaxExtent()==50);
    assert(vHead->GetNoSlices()==5);
// Check all nodes contain (correct) single volume
    for (G4int n=0;n<5;n++)
	{
	vProxy=vHead->GetSlice(n);
	assert(vProxy->IsNode());
	vNode=vProxy->GetNode();
	assert(vNode->GetNoContained()==1);
        assert(vNode->GetVolume(0)==n);
	}

// Open geometry freeing voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
// Check voxels destroyed
    assert(!myDaughter1Log.GetVoxelHeader());

    return true;
}
int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4SmartVoxelNodes());
    assert(testG4SmartVoxels());
    assert(testG4ReplicaVoxels());
    return 0;
}




