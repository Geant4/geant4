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
class BrachyMaterial;
class G4Sphere;
class G4Tubs;
class G4LogicalVolume;
class G4VisAttributes;

class BrachyDetectorConstructionLeipzig
{
public:
  BrachyDetectorConstructionLeipzig();
  ~BrachyDetectorConstructionLeipzig(); 
  void ConstructLeipzig(G4VPhysicalVolume*);
  void CleanLeipzigApplicator();

private:
  
  BrachyMaterial* pMaterial; 
  
  G4Tubs* capsule;
  G4Sphere* capsuleTip;
  G4Tubs* iridiumCore;   
  G4Tubs* applicator1;
  G4Tubs* applicator2;
  
  G4LogicalVolume* capsuleLog;
  G4LogicalVolume* capsuleTipLog;
  G4LogicalVolume* iridiumCoreLog;
  G4LogicalVolume* applicator1Log;
  G4LogicalVolume* applicator2Log;
 
  G4VPhysicalVolume* capsulePhys;
  G4VPhysicalVolume* capsuleTipPhys;
  G4VPhysicalVolume* iridiumCorePhys;
  G4VPhysicalVolume* applicator1Phys;
  G4VPhysicalVolume* applicator2Phys;

  G4VisAttributes* simpleCapsuleVisAtt;
  G4VisAttributes* simpleCapsuleTipVisAtt;
  G4VisAttributes* applicatorVisAtt;
};

#endif
