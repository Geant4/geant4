// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildBoxWorld.hh,v 1.2 1999-12-15 14:49:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef BUILDBOXWORLD_HH
#define BUILDBOXWORLD_HH

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"

G4VPhysicalVolume* BuildBoxWorld()
{
    G4double z = 13.;
    G4double a = 26.98*g/mole;
    G4double density = 2.7*g/cm3;

    G4Box *myWorldBox= new G4Box ("WBox",10000,10000,10000);
    G4Box *myTargetBox = new G4Box ("TBox",100,200,400);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						    "WLog",0,0,0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetBox,0,
						     "YLog",0,0,0);
    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(),
						  "TPhys",
						  myTargetLog,
						  myWorldPhys,false,0);
    
    return myWorldPhys;
}
#endif
