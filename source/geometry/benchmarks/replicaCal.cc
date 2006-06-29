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
// $Id: replicaCal.cc,v 1.8 2006-06-29 18:15:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Replicated geometry performance test - ~calorimeter with end caps. Main
// calorimeter block is replicated (divided) in r & z.
// Toggle optimisation with 0/1 on command line. Default - optimisation ON

#include "G4ios.hh"
#include <stdlib.h>
#include <cmath>

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

int main()
{
  G4ThreeVector pos(0,0,-10*x0+0.01*cm);
  G4double phi=pi*0.5,theta=15*deg;
  G4ThreeVector dir(std::cos(phi)*std::sin(theta),std::sin(phi)*std::sin(theta),std::cos(theta));
  G4VPhysicalVolume *myTopNode;
  G4Timer timer;

  G4cout << "***  Replication Performance Test - P.Kent 26.09.96  ***" << G4endl << G4endl;

  myTopNode=BuildReplicaCal(0);	// Build the geometry
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
  
  G4cout << G4endl << "Shooting from " << pos << " along " << dir << G4endl;
  ShootVerbose(myTopNode,pos,dir);
  Shoot(numShoot,myTopNode,pos,dir);

  G4GeometryManager::GetInstance()->OpenGeometry();
  return EXIT_SUCCESS;
}

