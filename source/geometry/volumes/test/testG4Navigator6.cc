// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Navigator6.cc,v 1.2 1999-11-24 21:13:54 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//   Create a tubular "calorimeter". Shoot from origin along x & y axes
//   printing location, steps & safeties. Locate 50^3 points within
//   calorimeter to check performance of point location logic.
//
//   Arguments: Put `1' or `0' to toggle voxel optimisation on/off
//              [Default: ON]
//
//   Define G4GEOMETRY_VERBOSE for dump of topmost voxels


#include <assert.h>
#include "G4ios.hh"
#include <stdlib.h>

// Global defs
#include "globals.hh"

#include "G4Timer.hh"
#include "ApproxEqual.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Build tubular calorimeter:
// An array of interlocking complete tubes, inside a box
//
// Controlled by:
const G4double kTubeHalfHeight = 10;
const G4double kTubeRadius = 5;
const G4double kTubeNoRow = 10;
const G4double kTubeNoColumn = 11; // Should be odd for symmetrical array

const G4double kBoxDx=kTubeNoRow*kTubeRadius;
const G4double yDelta=2.0*kTubeRadius*sin(M_PI/3.0);
const G4double kBoxDy=(kTubeNoColumn-1)*yDelta*0.5+kTubeRadius;
const G4double kBoxDz=kTubeHalfHeight;

G4VPhysicalVolume* BuildGeometry()
{
    G4double bigXStart=-(kTubeNoRow-1)*kTubeRadius;
    G4double smallXStart=bigXStart+kTubeRadius;

    G4double bigYStart=-(kTubeNoColumn-1)*yDelta*0.5;
    G4double smallYStart=bigYStart+yDelta;


    G4int row,column;

    G4Box *calBox = new G4Box ("Cal Box",kBoxDx,kBoxDy,kBoxDz);
    G4Tubs *calTube = new G4Tubs("Cal Tube",0,kTubeRadius,
				 kTubeHalfHeight,0,360);
    

    G4LogicalVolume *myDetectorLog=new G4LogicalVolume(calBox,0,
						       "World",
						       0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *myDetectorPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
						    "World",
						    myDetectorLog,0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume *calTubLog=new G4LogicalVolume(calTube,0,
						   "Cal Crystal",
						   0,0,0);

    G4String tname("Target");
    G4int copyNo=0;
    for (column=0;column<kTubeNoColumn;column+=2)
	{
	    for (row=0;row<kTubeNoRow;row++)
		{
		    
		    G4PVPlacement *calPhys=new G4PVPlacement(
			0,G4ThreeVector(bigXStart+row*kTubeRadius*2.0,bigYStart+column*yDelta,0),
			tname,calTubLog,
			myDetectorPhys,false,copyNo++);
		}
	}


    for (column=0;column<kTubeNoColumn-1;column+=2)
	{
	    for (row=0;row<kTubeNoRow-1;row++)
		{
		    
		    G4PVPlacement *calPhys=new G4PVPlacement(
			0,G4ThreeVector(smallXStart+row*kTubeRadius*2.0,smallYStart+column*yDelta),
			tname,calTubLog,
			myDetectorPhys,false,copyNo++);
		}
	}

    return myDetectorPhys;
}


G4bool printShoot(G4VPhysicalVolume *pTopNode,
		  const G4ThreeVector& pLoc,
		  const G4ThreeVector& pVec)
{
    G4double Step=0,safety=0;
    const G4double physStep=kInfinity;
    G4VPhysicalVolume *located=0;
    G4Navigator myNav;
    myNav.SetWorldVolume(pTopNode);

    G4ThreeVector partLoc(pLoc);
    G4cout << "Shooting from " << pLoc << " along " << pVec << endl;
    located=myNav.LocateGlobalPointAndSetup(partLoc);
    while (located)
	{
	    Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
	    G4cout << "Physical Location=" << located->GetName()
		 << " #" << located->GetCopyNo() << endl
		 << "   Step=" << Step << "  Safety=" << safety
		 << "  ---->" << endl;

	    partLoc+=Step*pVec;
	    myNav.SetGeometricallyLimitedStep();
	    located=myNav.LocateGlobalPointAndSetup(partLoc);
	}
    return true;
}

G4bool runLocate(G4VPhysicalVolume *pTopNode)
{
    const G4int numLocPerAxis=50;
    const G4double dxStep=kBoxDx*2.0/numLocPerAxis;
    const G4double dyStep=kBoxDy*2.0/numLocPerAxis;
    const G4double dzStep=kBoxDz*2.0/numLocPerAxis;

    G4Navigator myNav;
    G4ThreeVector worldPoint;
    myNav.SetWorldVolume(pTopNode);
    for (G4double x=-kBoxDx;x<kBoxDx;x+=dxStep)
	{
	    for (G4double y=-kBoxDy;y<kBoxDy;y+=dyStep)
		{
		    for (G4double z=-kBoxDz;z<kBoxDz;z+=dzStep)
			{
			    worldPoint=G4ThreeVector(x,y,z);
			    myNav.LocateGlobalPointAndSetup(worldPoint,0,false);
			}
		}
	    
	}
    return true;
}

G4bool runAll(G4VPhysicalVolume *pTopNode)
{
    G4cout << "Locating..." << endl;
    runLocate(pTopNode);
    G4cout << "Done" << endl;
    return true;
}

int main(int argc, char* argv[])
{
    G4bool optimise;
    G4Timer timer;
    if (argc==1)
	{
	    optimise=true;
	}
    else if (argc==2)
	{
	    G4String opt(argv[1]);

	    if (opt=="0")
		{
		    optimise=false;
		}
	    else if (opt=="1")
		{
		    optimise=true;
		}
	    else
		{
		    G4cout << "Unknown args" << endl;
		    return EXIT_FAILURE;
		}
	}
    else
	{
	    G4cout << "Unknown args" << endl;
	    return EXIT_FAILURE;
	}

    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry

    G4GeometryManager::GetInstance()->OpenGeometry();

    timer.Start();
    G4GeometryManager::GetInstance()->CloseGeometry(optimise);
    timer.Stop();

//#ifdef G4GEOMETRY_VERBOSE
//    G4cout << *(G4LogicalVolumeStore::GetInstance()->at(0)->GetVoxelHeader());
//#endif
    if (optimise)
	{
	    G4cout << "Built voxels ";
	}
    else
	{
	    G4cout << "No voxels ";
	}
    G4cout << timer << endl;

    printShoot(myTopNode,
 	       G4ThreeVector(-kBoxDx,0,0),
 	       G4ThreeVector(1,0,0));
     printShoot(myTopNode,
 	       G4ThreeVector(0,0,0),
	       G4ThreeVector(1/sqrt(2.),1/sqrt(2.),0));
    timer.Start();
    runAll(myTopNode);
    timer.Stop();
    G4cout << timer << endl;

    G4GeometryManager::GetInstance()->OpenGeometry();
    return EXIT_SUCCESS;
}





