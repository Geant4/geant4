// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: shooter.cc,v 1.1 1999-01-08 16:31:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// shooter - perform test shots through simple box world.
//
// World consisting of single box positioned inside large box world
// For comparison against F77 case. No voxels involved

#include "G4ios.hh"
#include <stdlib.h>

#include "G4GeometryManager.hh"
#include "BuildBoxWorld.hh"
#include "Shoot.hh"

const G4int numShoot = 10000;

int main()
{
    G4ThreeVector origin(0,0,0),pMX(-500,0,0);
    G4ThreeVector vx(1,0,0);
    G4ThreeVector vy(0,1,0);
    G4VPhysicalVolume *myTopNode;

    G4cout << "Simple BoxWorld Performance Test - P.Kent 21.08.95" << endl;
#ifndef NDEBUG
    G4cout << "WARNING: *** ASSERTs are compiled IN ***" << endl;
#endif

    myTopNode=BuildBoxWorld();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry();

    G4cout << endl << "Shooting from " << origin << " along " << vx << endl;
    ShootVerbose(myTopNode,origin,vx);
    Shoot(numShoot,myTopNode,origin,vx);

    G4cout << endl << "Shooting from " << origin << " along " << vy << endl;
    ShootVerbose(myTopNode,origin,vy);
    Shoot(numShoot,myTopNode,origin,vy);

    G4cout << endl << "Shooting from " << pMX << " along " << vx << endl;
    ShootVerbose(myTopNode,pMX,vx);
    Shoot(numShoot,myTopNode,pMX,vx);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return EXIT_SUCCESS;
}

