// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: replicaCal.cc,v 1.1 1999-01-08 16:31:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Replicated geometry performance test - ~calorimeter with end caps. Main
// calorimeter block is replicated (divided) in r & z.
// Toggle optimisation with 0/1 on command line. Default - optimisation ON

#include "G4ios.hh"
#include <stdlib.h>
#include <math.h>

#include "G4GeometryManager.hh"
#include "Shoot.hh"
#include "G4Timer.hh"

#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"

const G4int numShoot = 10000;
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

int main()
{
  G4ThreeVector pos(0,0,-10*x0+0.01*cm);
  G4double phi=M_PI*0.5,theta=15*deg;
  G4ThreeVector dir(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
  G4VPhysicalVolume *myTopNode;
  G4Timer timer;

  G4cout << "  Replication Performance Test - P.Kent 26.09.96" << endl;
#ifndef NDEBUG
  G4cout << "WARNING: *** ASSERTs are compiled IN ***" << endl;
#endif

  myTopNode=BuildReplicaCal(NULL);	// Build the geometry
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
  
  G4cout << endl << "Shooting from " << pos << " along " << dir << endl;
  ShootVerbose(myTopNode,pos,dir);
  Shoot(numShoot,myTopNode,pos,dir);

  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}

