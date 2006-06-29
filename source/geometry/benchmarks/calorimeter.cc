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
// $Id: calorimeter.cc,v 1.8 2006-06-29 18:15:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// calorimeter
//
// World consisting of CSoC '95 tube calorimeter
// For comparison against F77 case. Toggle optimisation with 0/1
// on command line. Default - optimisation ON

#include "G4ios.hh"
#include <stdlib.h>
#include <cmath>

#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
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
    G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);
    G4VPhysicalVolume *myTopNode;
    G4Timer timer;

    G4cout << "***  Calorimeter Performance Test - P.Kent 21.08.95  ***" << G4endl
	   << ">>>  Using calorimeter from CERN School of Computing 1995" << G4endl << G4endl;

    myTopNode=BuildCalorimeter();	// Build the geometry
    timer.Start();
    G4GeometryManager::GetInstance()->CloseGeometry(optimise);
    timer.Stop();

    if (optimise)
	{
	    G4cout << "Optimisation ON" << G4endl;
	}
    else
	{
	    G4cout << "Optimisation OFF" << G4endl;
	}
    G4cout << "Closing geometry took: " << timer << G4endl;

#ifdef G4GEOMETRY_VERBOSE
    G4cout << G4endl << "Voxels for:" 
	   << (*G4LogicalVolumeStore::GetInstance())[1]->GetName() << G4endl
	   << *(*G4LogicalVolumeStore::GetInstance())[1]->GetVoxelHeader()
	   << G4endl;
#endif

    G4SmartVoxelHeader* top=(*G4LogicalVolumeStore::GetInstance())[1]
	                    ->GetVoxelHeader();
    if (top)
	{
	    G4cout << "1st level voxels along ";
	    switch (top->GetAxis())
		{
		case kXAxis:
		    G4cout << "X" << G4endl;
		    break;
		case kYAxis:
		    G4cout << "Y" << G4endl;
		    break;
		case kZAxis:
		    G4cout << "Z" << G4endl;
		    break;
		default:
		    break;    
		}
	}

    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    ShootVerbose(myTopNode,origin,vx);
    Shoot(numShoot,myTopNode,origin,vx);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    ShootVerbose(myTopNode,origin,vy);
    Shoot(numShoot,myTopNode,origin,vy);

    G4cout << G4endl << "Shooting from " << origin << " along " << vxy << G4endl;
    ShootVerbose(myTopNode,origin,vxy);
    Shoot(numShoot,myTopNode,origin,vxy);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return EXIT_SUCCESS;
}

