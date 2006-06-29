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
// $Id: BuildCalorimeter.hh,v 1.6 2006-06-29 18:15:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef BUILDCALORIMETER_HH
#define BUILDCALORIMETER_HH

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

G4VPhysicalVolume* BuildCalorimeter()
{
    G4double offset=22.5,xTlate,yTlate;
    G4int i,j,copyNo;

    G4Box *myWorldBox= new G4Box ("WBox",10000.*mm,10000.*mm,10000.*mm);
    G4Box *myCalBox = new G4Box ("CBox",1500.*mm,1500.*mm,1000.*mm);
    G4Tubs *myTargetTube = new G4Tubs ("TTube",0.,22.5*mm,1000.*mm,0.,360.*deg);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						   "WLog",0,0,0);
    G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,0,
						  "CLog",0,0,0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,0,
						     "TLog",0,0,0);

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
    for (j=1;j<=25;j++)
      {
	yTlate=-1000.0-40.0+j*80.0;
	
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0-offset;
	    // G4PVPlacement *myTargetPhys=
	    new G4PVPlacement(0,G4ThreeVector(xTlate*mm,yTlate*mm,0),
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
    for (j=1;j<=26;j++)
      {
	yTlate=-1000.0-80.0+j*80.0;
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate=-1000.0-20.0+i*45.0;
	    // G4PVPlacement *myTargetPhys=
            new G4PVPlacement(0,G4ThreeVector(xTlate*mm,yTlate*mm,0),
					      tName2,
					      myTargetLog,
					      myCalPhys,
					      false,
					      copyNo);
	  }
      }

    return myWorldPhys;
}
#endif

