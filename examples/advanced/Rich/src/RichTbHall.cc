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
// Rich advanced example for Geant4
// RichTbHall.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include "RichTbDetectorConstruction.hh"
#include "RichTbHall.hh"
#include "RichTbMaterial.hh"
#include "RichTbGeometryParameters.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

RichTbHall::RichTbHall() { ; }
RichTbHall::RichTbHall(RichTbMaterial* RMaterial) {

  G4Box * RichTbHallBox
    = new G4Box("World",ExpHallHalfX,ExpHallHalfY,ExpHallHalfZ);
  G4LogicalVolume * RichTbHallLog
    = new G4LogicalVolume(RichTbHallBox,RMaterial->getAir(),"World",0,0,0);
  G4VPhysicalVolume * RichTbHallPhys
    = new G4PVPlacement(0,G4ThreeVector(),"World",RichTbHallLog,0,false,0);

  RichTbHallLVol  = RichTbHallLog;
  RichTbHallPVol  = RichTbHallPhys;
 


 }
RichTbHall::~RichTbHall() {; }
