// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: shooter.cc,v 1.3 2000-12-12 08:23:08 medernac Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// shooter - perform test shots through simple box world.
//
// World consisting of single box positioned inside large box world
// For comparison against F77 case. No voxels involved

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
  G4PVPlacement *leakPhys=new G4PVPlacement(0,G4ThreeVector(0,0,-zc),
					    "leakPhys",
					    leakLog,
					    ecalPhys,false,0);
  G4PVPlacement *leakPhys2=new G4PVPlacement(0,G4ThreeVector(0,0,zc),
					    "leakPhys",
					    leakLog,
					    ecalPhys,false,1);
  G4PVPlacement *latrPhys=new G4PVPlacement(0,G4ThreeVector(),
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
  G4PVReplica *blocPhysRZ=new G4PVReplica("blocRepRZ",
					  blocRZLog,blocPhysR,
					  kZAxis,nl,2*z1/nl);

  return ecalPhys;
}


/*
   to change accuracy :
	fieldMgr->GetChordFinder()->SetDeltaChord( G4double newValue);
 */
/**
   Speed tests and benchmarking
   
   Options:
   - Number of events
   - Geometry type ("box": simple box,"caloloop": calorimeter with a creation loop,
   "calorep": calorimeter done by replication)
   - Magnetic Field  [value of the field, accuracy] (default OFF)
**/
int main(int argc,char *argv[])
{
  G4ThreeVector origin(0,0,0),pMX(-500,0,0);
  G4ThreeVector vx(1,0,0);
  G4ThreeVector vy(0,1,0);
  G4VPhysicalVolume *myTopNode;

  G4double Field ;
  G4double DeltaChord ;

  enum GeomType { BOX, CALOLOOP, CALOREP} GeomType;
  G4int NOMAGFIELD ;

  int i;

  /* Default Value */
  
  GeomType = BOX ;
  NOMAGFIELD = 1 ;
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
      
      NOMAGFIELD = 0;
    }
  }
  
  
  G4cout << "Magnetic Navigation Performace tester -- E.Medernach 30/10/00" << G4endl;
#ifndef NDEBUG
  G4cout << "WARNING: *** ASSERTs are compiled IN ***" << G4endl << G4endl ;
#endif

  // Build the geometry
  
  switch (GeomType) {

  case BOX:
    G4cout << " Box only " << G4endl ;
    myTopNode=BuildBoxWorld();
    break;

  case CALOLOOP:
    G4cout << " Calorimeter with loop " << G4endl ;
    myTopNode=BuildCalorimeter();
    break;

  case CALOREP:
    G4cout << "  Calorimeter with replica " << G4endl ;
    myTopNode=BuildReplicaCal(NULL);
    break;
    
  }
  
   
  G4GeometryManager::GetInstance()->CloseGeometry(true);

  if (NOMAGFIELD) {

    G4cout << "No Magnetic Field test" << G4endl ;
    
    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    ShootVerbose(myTopNode,origin,vx);
    Shoot(numShoot,myTopNode,origin,vx);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    ShootVerbose(myTopNode,origin,vy);
    Shoot(numShoot,myTopNode,origin,vy);

    G4cout << G4endl << "Shooting from " << pMX << " along " << vx << G4endl;
    ShootVerbose(myTopNode,pMX,vx);
    Shoot(numShoot,myTopNode,pMX,vx);

  } else {

    //    Field = 0.1 * tesla ;
    DeltaChord = 1.0e-2 * mm ;
    
    
    G4cout << "Magnetic Field test with " << Field/tesla << " Tesla" << G4endl ;

    G4cout << G4endl << "Shooting from " << origin << " along " << vx << G4endl;
    //    ShootVerbose(myTopNode,origin,vx);
    MagneticShoot(numShoot,myTopNode,origin,vx,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << origin << " along " << vy << G4endl;
    //    ShootVerbose(myTopNode,origin,vy);
    MagneticShoot(numShoot,myTopNode,origin,vy,Field,DeltaChord);

    G4cout << G4endl << "Shooting from " << pMX << " along " << vx << G4endl;
    //    ShootVerbose(myTopNode,pMX,vx);
    MagneticShoot(numShoot,myTopNode,pMX,vx,Field,DeltaChord);

  }
    
    
    
  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}



