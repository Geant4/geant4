// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4NavigationHistory.cc,v 1.2 1999-11-24 21:13:03 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// testG4NavigationHistory
//             Ensure asserts are compiled in
//

#include <assert.h>

// Global defs
#include "globals.hh"

// Tested entities
#include "G4NavigationHistory.hh"

// Required entities
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Combined test of logical and physical volumes
// 2 small cubes are positioned inside a larger cuboid
//
// Check navigation links, `Store' entries
G4bool testG4NavigationHistory()
{
    G4Box myBigBox("cuboid",25,20,20);
    G4Box myBox("cube",10,10,10);
    G4ThreeVector vorigin,vmx(-15,0,0),vx(15,0,0);

    G4LogicalVolume detectorLog(&myBigBox,0,
				"World",0,0,0);
				// Logical with no mag field,
                                // sensitive detector or user limits
    
    G4PVPlacement detectorPhys(0,vorigin,
			       "World",&detectorLog,
			       0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume myDaughterLog(&myBox,0,"block",0,0,0);
    G4PVPlacement offXPhys(0,vx,
			   "Target 1",&myDaughterLog,
			  &detectorPhys,false,0);
    G4PVPlacement offMXPhys(0,vmx,
			    "Target 2",&myDaughterLog,
			    &detectorPhys,false,0);

// 
// BEGIN test of G4NaviagtionHistory
//
    G4NavigationHistory nHist;

    nHist.SetFirstEntry(&detectorPhys);
    assert(nHist.GetDepth()==0);
    assert(nHist.GetMaxDepth()>0);
    assert(nHist.GetTransform(0)==G4AffineTransform(vorigin));
    assert(nHist.GetReplicaNo(0)==-1);
    assert(nHist.GetVolume(0)==&detectorPhys);

    offXPhys.Setup(&detectorPhys);
    nHist.NewLevel(&offXPhys);
    assert(nHist.GetDepth()==1);
    assert(nHist.GetTransform(1)==G4AffineTransform(vmx));
    assert(nHist.GetReplicaNo(1)==-1);
    assert(nHist.GetVolume(0)==&detectorPhys);
    assert(nHist.GetVolume(1)==&offXPhys);

    nHist.BackLevel();
    assert(nHist.GetDepth()==0);
    assert(nHist.GetTransform(0)==G4AffineTransform(vorigin));
    assert(nHist.GetReplicaNo(0)==-1);
    assert(nHist.GetVolume(0)==&detectorPhys);

    offMXPhys.Setup(&detectorPhys);
    nHist.NewLevel(&offMXPhys);
    assert(nHist.GetDepth()==1);
    assert(nHist.GetTransform(1)==G4AffineTransform(vx));
    assert(nHist.GetReplicaNo(1)==-1);
    assert(nHist.GetVolume(0)==&detectorPhys);
    assert(nHist.GetVolume(1)==&offMXPhys);
    
//
// END test of G4NaviagtionHistory
//
    return true;
}

G4int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4NavigationHistory());
    return 0;
}





