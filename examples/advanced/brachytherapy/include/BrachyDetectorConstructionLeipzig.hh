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
//
//
//    ************************************************
//    *                                              *
//    *      BrachyDetectorConstructionLeipzig.hh    *
//    *                                              *
//    ************************************************
//
// Management and modelling of Leipzig applicator 
//
// Code by S. Guatelli
//
#ifndef BrachyDetectorConstructionLeipzig_H
#define BrachyDetectorConstructionLeipzig_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4Sphere;
class G4Tubs;
class G4LogicalVolume;
class G4VisAttributes;

class BrachyDetectorConstructionLeipzig
{
public:
  BrachyDetectorConstructionLeipzig();
  ~BrachyDetectorConstructionLeipzig()=default; 
  void ConstructLeipzig(G4VPhysicalVolume*);
  void CleanLeipzigApplicator();

private:
  G4Tubs* fCapsule;
  G4Sphere* fCapsuleTip;
  G4Tubs* fIridiumCore;   
  G4Tubs* fApplicator1;
  G4Tubs* fApplicator2;
  
  G4LogicalVolume* fCapsuleLog;
  G4LogicalVolume* fCapsuleTipLog;
  G4LogicalVolume* fIridiumCoreLog;
  G4LogicalVolume* fApplicator1Log;
  G4LogicalVolume* fApplicator2Log;
 
  G4VPhysicalVolume* fCapsulePhys;
  G4VPhysicalVolume* fCapsuleTipPhys;
  G4VPhysicalVolume* fIridiumCorePhys;
  G4VPhysicalVolume* fApplicator1Phys;
  G4VPhysicalVolume* fApplicator2Phys;

  G4VisAttributes* fSimpleCapsuleVisAtt;
  G4VisAttributes* fSimpleCapsuleTipVisAtt;
  G4VisAttributes* fApplicatorVisAtt;
};
#endif
