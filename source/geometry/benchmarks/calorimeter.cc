// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: calorimeter.cc,v 1.1 1999-01-08 16:31:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// calorimeter
//
// World consisting of CSoC '95 tube calorimeter
// For comparison against F77 case. Toggle optimisation with 0/1
// on command line. Default - optimisation ON

#include "G4ios.hh"
#include <stdlib.h>
#include <math.h>

#include "G4GeometryManager.hh"
#include "BuildCalorimeter.hh"
#include "Shoot.hh"
#include "G4Timer.hh"
#include "globals.hh"
#include "geomdefs.hh"

const G4int numShoot = 10000;
const G4bool optimise= true;

int main()
{
    G4ThreeVector origin(0,0,0),pMX(-500,0,0);
    G4ThreeVector vx(1,0,0);
    G4ThreeVector vy(0,1,0);
    G4ThreeVector vxy(1/sqrt(2.0),1/sqrt(2.0),0);
    G4VPhysicalVolume *myTopNode;
    G4Timer timer;

    G4cout << "  Calorimeter Performance Test - P.Kent 21.08.95" << endl
	 << "Using calorimeter from CERN School of Computing '95" << endl;
#ifndef NDEBUG
    G4cout << "WARNING: *** ASSERTs are compiled IN ***" << endl;
#endif

    myTopNode=BuildCalorimeter();	// Build the geometry
    timer.Start();
    G4GeometryManager::GetInstance()->CloseGeometry(optimise);
    timer.Stop();

    if (optimise)
	{
	    G4cout << "Optimisation ON";
	}
    else
	{
	    G4cout << "Optimisation OFF";
	}
    G4cout << " Geometry close took " << timer << endl;

#ifdef G4GEOMETRY_VERBOSE
    G4cout << endl << "Voxels for:" 
	 << G4LogicalVolumeStore::GetInstance()->at(1)->GetName() << endl
	<< *G4LogicalVolumeStore::GetInstance()->at(1)->GetVoxelHeader()
	<< endl;
#endif
    G4SmartVoxelHeader* top=G4LogicalVolumeStore::GetInstance()
	                   ->at(1)->GetVoxelHeader();
    if (top)
	{
	    G4cout << "1st level voxels along ";
	    switch (top->GetAxis())
		{
		case kXAxis:
		    G4cout << "X" << endl;
		    break;
		case kYAxis:
		    G4cout << "Y" << endl;
		    break;
		case kZAxis:
		    G4cout << "Z" << endl;
		    break;
		    
		}
	}

    G4cout << endl << "Shooting from " << origin << " along " << vx << endl;
    ShootVerbose(myTopNode,origin,vx);
    Shoot(numShoot,myTopNode,origin,vx);

    G4cout << endl << "Shooting from " << origin << " along " << vy << endl;
    ShootVerbose(myTopNode,origin,vy);
    Shoot(numShoot,myTopNode,origin,vy);

    G4cout << endl << "Shooting from " << origin << " along " << vxy << endl;
    ShootVerbose(myTopNode,origin,vxy);
    Shoot(numShoot,myTopNode,origin,vxy);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return EXIT_SUCCESS;
}

