//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: shooter.cc,v 1.5 2002-01-09 16:17:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// shooter - perform test shots.
//

/*
  gmake CXXFLAGS="-g -pg -a -lc_p " CPPVERBOSE=1 for the library

  gmake CXXFLAGS="-g -pg -a -lc_p -static " CPPVERBOSE=1 G4TARGET=shooter
  shooter
  (line by line profiling in detail)
  gprof -i -p -q -x -A -l `which shooter` > profile.shooter
  (normal profiling)
  gprof -p `which shooter` > profile.shooter
*/

#include "G4ios.hh"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "G4GeometryManager.hh"
#include "BuildBoxWorld.hh"
#include "BuildCalorimeter.hh"
#include "Shoot.hh"

#include "G4GeometryManager.hh"
#include "Shoot.hh"
#include "G4Timer.hh"

#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"


G4int numShoot ; //1000000;
  

const G4bool optimise= true;
const G4double x0=1.12343*cm;

G4VPhysicalVolume* BuildReplicaCal(G4Material* Air)
{
  G4double xr=x0/4;
  G4int nr=20;
  G4int nl=20;

  G4double r1=20*xr,r2=21*xr,z1=10*x0,z2=11*x0;
  // Container
  G4Tubs *ecalTube=new G4Tubs("ECAL",0,r2,z2,0,360*deg);
  // End cap
  G4Tubs *leakTube=new G4Tubs("LEAK",0,r2,0.5*x0,0,360*deg);
  // Wrapper
  G4Tubs *latrTube=new G4Tubs("LATR",r1,r2,z1,0,360*deg);
  // Main calorimeter block
  G4Tubs *blocTube=new G4Tubs("BLOC",0,r1,z1,0,360*deg);
  G4Tubs *blocTubeR=new G4Tubs("BLOCR",0,r1/nr,z1,0,360*deg);
  G4Tubs *blocTubeRZ=new G4Tubs("BLOCRZ",0,r1/nr,z1/nl,0,360*deg);

  G4double zc=0.5*(z1+z2);

  G4LogicalVolume *ecalLog=new G4LogicalVolume(ecalTube,Air,
					       "ecalLog",0,0,0);
  G4LogicalVolume *leakLog=new G4LogicalVolume(leakTube,Air,
					       "leakLog",0,0,0);
  G4LogicalVolume *latrLog=new G4LogicalVolume(latrTube,Air,
					       "latrLog",0,0,0);
  G4LogicalVolume *blocLog=new G4LogicalVolume(blocTube,Air,
					       "blocLog",0,0,0);
  G4LogicalVolume *blocRLog=new G4LogicalVolume(blocTubeR,Air,
					       "blocRLog",0,0,0);
  G4LogicalVolume *blocRZLog=new G4LogicalVolume(blocTubeRZ,Air,
						 "blocRZLog",0,0,0);


  G4PVPlacement *ecalPhys=new G4PVPlacement(0,G4ThreeVector(),
					    "ecalPhys",
					    ecalLog,
					    0,false,0);
  // Position end caps, wrapper and calorimeter bloc within ecal
  // G4PVPlacement *leakPhys=
  new G4PVPlacement(0,G4ThreeVector(0,0,-zc),
		    "leakPhys",
		    leakLog,
		    ecalPhys,false,0);
  // G4PVPlacement *leakPhys2=
  new G4PVPlacement(0,G4ThreeVector(0,0,zc),
		    "leakPhys",
		    leakLog,
		    ecalPhys,false,1);
  // G4PVPlacement *latrPhys=
  new G4PVPlacement(0,G4ThreeVector(),
		    "latrPhys",
		    latrLog,
		    ecalPhys,false,0);
  G4PVPlacement *blocPhys=new G4PVPlacement(0,G4ThreeVector(),
					    "blocPhys",
					    blocLog,
					    ecalPhys,false,0);

  // Create replicas
  G4PVReplica *blocPhysR=new G4PVReplica("blocRepR",
					 blocRLog,blocPhys,
					 kRho,nr,r1/nr,0);
  // G4PVReplica *blocPhysRZ=
  new G4PVReplica("blocRepRZ",
		  blocRZLog,blocPhysR,
		  kZAxis,nl,2*z1/nl);

  return ecalPhys;
}


/*
   to change accuracy :
	fieldMgr->GetChordFinder()->SetDeltaChord( G4double newValue);
*/

int main(int argc,char *argv[])
{
  G4ThreeVector origin(0,0,0),pMX(-500,0,0);
  G4ThreeVector vx(1,0,0);
  G4ThreeVector vy(0,1,0);
  G4VPhysicalVolume *myTopNode=0;

  G4double Field = 0.*tesla;

  enum GeomType { BOX, CALOLOOP, CALOREP} GeomType;
  G4int i;
  G4bool useField = false;

  /* Default Value */
  
  GeomType = BOX ;
  numShoot = 1000000 ;
  
  /* Command line parsing */
  
  for (i=1;i<argc;i++) {
    if ((i < (argc-1)) && (strcmp(argv[i],"-event") == 0)) {
      sscanf(argv[i+1],"%d",&numShoot);      
    }
    if ((i < (argc-1)) && (strcmp(argv[i],"-geom") == 0)) {
      if (strcmp(argv[i+1],"box") == 0) {
	GeomType = BOX;
      } else    
      if (strcmp(argv[i+1],"caloloop") == 0) {
	GeomType = CALOLOOP;
      } else      
      if (strcmp(argv[i+1],"calorep") == 0) {
	GeomType = CALOREP;
      } else {
	G4cerr << argv[i+1] << " is not a known geometry (box,caloloop,calorep)" << G4endl ;
	exit (0);
	
      }
    }
    if ((i < (argc-1)) && (strcmp(argv[i],"-magn") == 0)) {
      sscanf(argv[i+1],"%lf",&Field);
      G4cout << " Mag Field = " << Field << G4endl ;
      
      Field = Field * tesla ;
      useField = true;
    }
  }
  
  
  G4cout << "***  Navigation Performance Tester - E.Medernach 30.10.00  ***" << G4endl
	 << ">>>  Based on original benchmark test by P.Kent" << G4endl << G4endl;

  G4cout << "Options (as arguments):" << G4endl
         << "-event <number_of_events>" << G4endl
         << "      number of events for the test. Default is 1000000" << G4endl
         << "-geom <geometry_type>" << G4endl
         << "      where <geometry_type> can be:" << G4endl
         << "      box      - simple box (default)" << G4endl
         << "      caloloop - calorimeter made by a loop of placements" << G4endl
         << "      calorep  - calorimeter made of replicas" << G4endl
         << "-magn <magnetic_field_value>" << G4endl
         << "      activates magnetic field (value in tesla units). Default is OFF" << G4endl << G4endl;

  // Build the geometry

  G4cout << "Geometry type:";  
  switch (GeomType) {

  case BOX:
    G4cout << " Box only." << G4endl ;
    myTopNode=BuildBoxWorld();
    break;

  case CALOLOOP:
    G4cout << " Calorimeter made of placements." << G4endl ;
    myTopNode=BuildCalorimeter();
    break;

  case CALOREP:
    G4cout << "  Calorimeter made of replicas." << G4endl ;
    myTopNode=BuildReplicaCal(0);
    break;
    
  }
  
  G4GeometryManager::GetInstance()->CloseGeometry(true);

  if (!useField)
  {
    G4cout << "--> Magnetic Field is disabled !" << G4endl ;
    
    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    ShootVerbose(myTopNode,origin,vx);
    Shoot(numShoot,myTopNode,origin,vx);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    ShootVerbose(myTopNode,origin,vy);
    Shoot(numShoot,myTopNode,origin,vy);

    G4cout << G4endl << "Shooting from " << pMX << " along " << vx << G4endl;
    ShootVerbose(myTopNode,pMX,vx);
    Shoot(numShoot,myTopNode,pMX,vx);
  }
  else
  {
    //    Field = 0.1 * tesla ;
    G4double DeltaChord = 1.0e-2 * mm ;
    
    G4cout << "--> Magnetic Field test with " << Field/tesla << " Tesla" << G4endl ;

    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    ShootVerbose(myTopNode,origin,vx);
    MagneticShoot(numShoot,myTopNode,origin,vx,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    ShootVerbose(myTopNode,origin,vy);
    MagneticShoot(numShoot,myTopNode,origin,vy,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << pMX << " along " << vx << G4endl;
    ShootVerbose(myTopNode,pMX,vx);
    MagneticShoot(numShoot,myTopNode,pMX,vx,Field,DeltaChord);
  }

  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}

