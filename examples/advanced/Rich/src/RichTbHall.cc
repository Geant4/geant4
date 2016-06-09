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
// Rich advanced example for Geant4
// RichTbHall.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
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
