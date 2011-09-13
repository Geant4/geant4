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
// $Id: testG4Navigator6.cc,v 1.5 2006-06-29 18:37:28 gunter Exp $
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

// Build wedge calorimeter:
// An array of tubs wedges, sharing the Z axis.
//
// Controlled by:
const G4double kTubeHalfHeight = 100.0 * mm;
const G4double kTubeRadius =   10.0 * mm ;
const G4int    kNumberOfWedges = 25;

const G4double kBoxDx=1.1 * kTubeRadius;
// const G4double yDelta=2.0*kTubeRadius*std::sin(pi/3.0);
const G4double kBoxDy=1.1 * kTubeRadius; 
const G4double kBoxDz=1.1 * kTubeHalfHeight;

G4RotationMatrix rotmWedge[kNumberOfWedges]; 

G4VPhysicalVolume* BuildGeometry()
{
    G4Box *calBox = new G4Box ("Cal Box",kBoxDx,kBoxDy,kBoxDz);
    G4Tubs *calTube = new G4Tubs("Wedge",
				 0,
				 kTubeRadius,
				 kTubeHalfHeight,
				 0,
				 (360.0*degree)/kNumberOfWedges 
				);
    G4cout << " Created Tubs solid with radius " << kTubeRadius
	   << " and halfHeight " << kTubeHalfHeight << G4endl;

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

    G4String tname("Wedge");
    G4int copyNo=0;
    G4ThreeVector  axisZ( 0.0, 0.0, 1.0 ); 
    G4int phiNum;

    for (phiNum=0;phiNum<kNumberOfWedges;phiNum++)
    {
        G4double angle= (phiNum+0.5) * 360.0 * degree / kNumberOfWedges; 
        rotmWedge[ phiNum ] = G4RotationMatrix( axisZ, angle ); 
       
//	   G4PVPlacement *calPhys=
        new G4PVPlacement(
			  &( rotmWedge[phiNum] ), 
			  G4ThreeVector( 0., 0., 0.),
			  tname,
			  calTubLog,
			  myDetectorPhys,
			  false,
			  copyNo++);
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
    MyNavigator myNav;
    myNav.SetWorldVolume(pTopNode);

    G4ThreeVector partLoc(pLoc);
    G4cout << G4endl;
    G4cout << "Shooting from " << pLoc << " along " << pVec << G4endl;
    located=myNav.LocateGlobalPointAndSetup(partLoc);
    unsigned int  stepNo=0; 
    while (located)
    {
       stepNo++; 
       Step=myNav.ComputeStep(partLoc, pVec, physStep, safety);
       G4cout << " Step# " << std::setw(4) << stepNo << "  ";
       G4cout << "Physical Location=" << std::setw(15) << located->GetName()
	      << " #" << std::setw(4) << located->GetCopyNo() 
	      << "   Step=" << std::setw(10) << Step 
	      << " Safety=" << std::setw(10) << safety;
       G4cout << "  ---->" 
	      << std::setw(14) << partLoc 
	      << G4endl;
       
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

    MyNavigator myNav;
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
    G4cout << "Locating..." << G4endl;
    runLocate(pTopNode);
    G4cout << "Done" << G4endl;
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
		    G4cout << "Unknown args" << G4endl;
		    return EXIT_FAILURE;
		}
	}
    else
	{
	    G4cout << "Unknown args" << G4endl;
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
    G4cout << timer << G4endl;

    printShoot(myTopNode,
 	       G4ThreeVector(0.,0.,0.),
 	       G4ThreeVector(0,0,1.0));
    printShoot(myTopNode,
 	       G4ThreeVector(0,0,0),
	       G4ThreeVector(1/std::sqrt(2.),1/std::sqrt(2.),0));
    timer.Start();
    runAll(myTopNode);
    timer.Stop();
    G4cout << timer << G4endl;

    G4GeometryManager::GetInstance()->OpenGeometry();
    return EXIT_SUCCESS;
}
