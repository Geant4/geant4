// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildBoxWorld.hh,v 1.3 2000-12-12 08:23:07 medernac Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef BUILDBOXWORLD_HH
#define BUILDBOXWORLD_HH

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Material.hh"

G4VPhysicalVolume* BuildBoxWorld()
{
    G4double z = 13.;
    G4double a = 26.98*g/mole;
    G4double density = 2.7*g/cm3;

    G4Box *myWorldBox= new G4Box ("WBox",10000*mm,10000*mm,10000*mm);
    G4Box *myTargetBox = new G4Box ("TBox",100*mm,200*mm,400*mm);
    //G4Cons *myTargetBox = new G4Cons ("TBox",100*mm,200*mm,400*mm,600*mm,200*mm,0,2*M_PI);
 
    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						    "WLog",0,0,0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetBox,0,
						     "YLog",0,0,0);
    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(10*mm,10*mm,10*mm),
						  "TPhys",
						  myTargetLog,
						  myWorldPhys,false,0);
    
    return myWorldPhys;
}
#endif
