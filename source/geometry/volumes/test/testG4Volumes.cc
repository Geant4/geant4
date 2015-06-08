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
// $Id$
//
// 
// testG4Volumes Placement and replica/parameterised volumes setup tests
//               Ensure asserts are compiled in

#include <assert.h>

#include "globals.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4Box.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Test of logical volume, where not related to physical volumes
G4bool testG4LogicalVolume()
{
    G4String myName("MySolenoid");
    G4String myName2("red");
    G4Box myBox("cube",10,10,10);
    G4Box myBigBox("cuboid",25,20,20);

    G4LogicalVolume myMotherVol(&myBigBox,0,myName,0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
// Check naming
    assert(myMotherVol.GetName()==myName);
    myMotherVol.SetName(myName2);
    assert(myMotherVol.GetName()==myName2);
// Check remaining initialisation
    assert(myMotherVol.GetNoDaughters()==0);
    assert(myMotherVol.GetSolid()==&myBigBox);
    assert(myMotherVol.GetMaterial()==0);
    assert(myMotherVol.GetFieldManager()==0);
    assert(myMotherVol.GetSensitiveDetector()==0);
    assert(myMotherVol.GetUserLimits()==0);
    assert(myMotherVol.GetVoxelHeader()==0);

// Check daughter access functionality
    G4LogicalVolume myDaughter(&myBox,0,"block",0,0,0);
    G4PVPlacement offsetx(0,G4ThreeVector(-15,0,0),
			  "Target 1",&myDaughter,
			  0,false,0);
				// Note: no mother pointer set
    G4PVPlacement offsety(0,G4ThreeVector(15,0,0),
			  "Target 2",&myDaughter,
			  0,false,0);
				// Note: no mother pointer set

    assert(myMotherVol.GetNoDaughters()==0);
    assert(myMotherVol.IsDaughter(&offsetx)==false);

    myMotherVol.AddDaughter(&offsetx);
    assert(myMotherVol.GetNoDaughters()==1);
    assert(myMotherVol.IsDaughter(&offsetx)==true);
    assert(myMotherVol.IsDaughter(&offsety)==false);
    assert(myMotherVol.GetDaughter(0)==&offsetx);
    myMotherVol.RemoveDaughter(&offsetx);
    assert(myMotherVol.GetNoDaughters()==0);
    assert(myMotherVol.IsDaughter(&offsetx)==false);
// myMotherVol left unchanged

    return true;
}

// Test of simple positioning, where not related to logical volumes
G4bool testG4PVPlacement()
{
    G4Box myBox("cube",10,10,10);
    G4ThreeVector vmx(-15,0,0);
    G4LogicalVolume myDaughter(&myBox,0,"block",0,0,0);
    G4PVPlacement offsetx(0,vmx,"Target 1",&myDaughter,
			  0,false,0);

    assert(offsetx.GetCopyNo()==0);
    assert(offsetx.GetTranslation()==vmx);
    assert(offsetx.GetRotation()==0);
    assert(offsetx.GetLogicalVolume()==&myDaughter);
    assert(offsetx.GetName()=="Target 1");
    assert(offsetx.IsMany()==false);
    assert(offsetx.IsReplicated()==false);

    assert(offsetx.GetTranslation()==vmx);
    assert(offsetx.GetRotation()==0);

    assert(offsetx.GetFrameRotation()==0);
    G4cout << "FrameTranslation=" << offsetx.GetFrameTranslation() 
         << " - Expected= "        << vmx << G4endl;
//    assert(offsetx.GetFrameTranslation()==vmx);

    assert(offsetx.GetObjectTranslation()==vmx);
    assert(offsetx.GetObjectRotation()->isIdentity());
    assert(offsetx.GetObjectRotationValue().isIdentity());

    return true;
}

// Test of simple replication, where not related to logical volumes
G4bool testG4PVReplica()
{
    G4Box myBigBox("Big Cube",1000,1000,1000);
    G4LogicalVolume myMotherVol(&myBigBox,0,"Big cube",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits

    G4Box myBox("cube",10,10,10);
    G4LogicalVolume myDaughter(&myBox,0,"block",0,0,0);
    G4PVReplica replicaX("Target 1",&myDaughter,&myMotherVol,
			kXAxis,2,10);

    assert(replicaX.GetLogicalVolume()==&myDaughter);
    assert(replicaX.GetName()=="Target 1");
    assert(replicaX.IsMany()==false);
    assert(replicaX.GetCopyNo()== (-1));
    assert(replicaX.IsReplicated()==true);

    EAxis axis;
    G4double width,offset;
    G4int n;
    G4bool consuming;
    replicaX.GetReplicationData(axis,n,width,offset,consuming);
    assert(axis==kXAxis&&n==2&&width==10&&offset==0&&consuming==true);

    return true;
}

// Combined test of logical and physical volumes
// 2 small cubes are positioned inside a larger cuboid
//
// Check navigation links, `Store' entries
G4bool testG4Volumes()
{
    G4Box myBigBox("cuboid",25,20,20);
    G4Box myBox("cube",10,10,10);

    G4LogicalVolume myDetectorLog(&myBigBox,0,
				  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement myDetectorPhys(0,G4ThreeVector(0,0,0),
				"World",&myDetectorLog,
				0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume myDaughterLog(&myBox,0,"block",0,0,0);
    G4PVPlacement offXPhys(0,G4ThreeVector(-15,0,0),
			  "Target 1",&myDaughterLog,
			  &myDetectorPhys,false,0);
    G4PVPlacement offMXPhys(0,G4ThreeVector(15,0,0),
			      "Target 2",&myDaughterLog,
			      &myDetectorPhys,false,0);

// Assert on navigation links
    assert(offXPhys.GetMotherLogical()==&myDetectorLog);
    assert(offMXPhys.GetMotherLogical()==&myDetectorLog);

    assert(myDetectorPhys.GetLogicalVolume()==&myDetectorLog);
    assert(myDetectorLog.GetNoDaughters()==2);
    assert(myDetectorLog.GetDaughter(0)==&offXPhys);
    assert(myDetectorLog.GetDaughter(1)==&offMXPhys);

// Check stores are ok
    assert(G4PhysicalVolumeStore::GetInstance()->size()==3);
    assert(G4LogicalVolumeStore::GetInstance()->size()==2);
    G4SolidStore* solidStore = G4SolidStore::GetInstance();
    G4bool exists = 0;
    std::vector<G4VSolid*>::const_iterator i;
    assert(solidStore->size()==2);
    for (i=solidStore->begin(); i!=solidStore->end(); ++i) {
      if (**i==myBox) {
	exists = 1;
	break;
      }
    }
    assert(exists);

// Check setup set mother ok
    assert(offXPhys.GetMotherLogical()==&myDetectorLog);
    assert(offXPhys.GetTranslation()==G4ThreeVector(-15,0,0));
    assert(offXPhys.GetRotation()==0);

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("main","000",FatalException,"FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4LogicalVolume());
    assert(testG4PVPlacement());
    assert(testG4PVReplica());
    assert(testG4Volumes());
    return 0;
}





