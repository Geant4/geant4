// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildCalorimeter.cc,v 1.2 1999-12-15 14:54:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "BuildCalorimeter.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

G4VPhysicalVolume* BuildCalorimeter()
{
    G4double offset=22.5,xTlate,yTlate;
    G4int i,j,copyNo;
    G4double z = 13.;
    G4double a = 26.98*g/mole;
    G4double density = 2.7*g/cm3;
    G4Material myMaterial("Aluminium", z, a, density);

    G4VisAttributes* red   = new G4VisAttributes (G4Colour(1,0,0));
    G4VisAttributes* green = new G4VisAttributes (G4Colour(0,1,0));
    G4VisAttributes* blue  = new G4VisAttributes (G4Colour(0,0,1));

    G4Box *myWorldBox= new G4Box ("WBox",10000,10000,10000);
    G4Box *myCalBox = new G4Box ("CBox",1500,1500,1000);
//    G4Tubs *myTargetTube = new G4Tubs ("TTube",0,22.5,1000,0,360);
    G4Tubs *myTargetTube = new G4Tubs ("TTube",0,22.5,100,0,360);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,&myMaterial,
						    "WLog",0,0,0);
    G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,&myMaterial,
						  "CLog",0,0,0);
    myCalLog -> SetVisAttributes (green);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,&myMaterial,
						     "TLog",0,0,0);
    myTargetLog -> SetVisAttributes (blue);

    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
					       myWorldPhys,false,0);

    G4String tName1("TPhys1");	// Allow all target physicals to share
    // same name (delayed copy)
    copyNo=0;
    //for (j=1;j<=25;j++)
    for (j=1;j<=3;j++)
      {
	yTlate=-1000.0-40.0+j*80.0;
	
	//for (i=1;i<=50;i++)
	for (i=1;i<=5;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0-offset;
	    G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0),
							  tName1,
							  myTargetLog,
							  myCalPhys,
							  false,
							  copyNo);
	  }
      }

    G4String tName2("TPhys2");	// Allow all target physicals to share
				// same name (delayed copy)
    copyNo=0;
    //for (j=1;j<=26;j++)
    for (j=1;j<=4;j++)
      {
	yTlate=-1000.0-80.0+j*80.0;
	//for (i=1;i<=50;i++)
	for (i=1;i<=5;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0;
	    G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0),
							  tName2,
							  myTargetLog,
							  myCalPhys,
							  false,
							  copyNo);
	  }
      }

    return myWorldPhys;
}
