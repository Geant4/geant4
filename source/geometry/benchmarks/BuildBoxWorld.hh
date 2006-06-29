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
// $Id: BuildBoxWorld.hh,v 1.7 2006-06-29 18:15:24 gunter Exp $
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
    G4Box *myWorldBox= new G4Box ("WBox",10000*mm,10000*mm,10000*mm);
    G4Box *myTargetBox = new G4Box ("TBox",100*mm,200*mm,400*mm);
    //G4Cons *myTargetBox = new G4Cons ("TBox",100*mm,200*mm,400*mm,600*mm,200*mm,0,2*pi);
 
    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,0,
						    "WLog",0,0,0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetBox,0,
						     "YLog",0,0,0);
    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    // G4PVPlacement *myTargetPhys=
    new G4PVPlacement(0,G4ThreeVector(10*mm,10*mm,10*mm),
		      "TPhys",
		      myTargetLog,
		      myWorldPhys,false,0);
    
    return myWorldPhys;
}
#endif
