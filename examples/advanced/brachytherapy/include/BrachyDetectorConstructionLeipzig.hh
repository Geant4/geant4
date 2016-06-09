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
//
// $Id: BrachyDetectorConstructionLeipzig.hh,v 1.3 2003/05/22 17:20:41 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
//    ************************************************
//    *                                              *
//    *      BrachyDetectorConstructionLeipzig.hh    *
//    *                                              *
//    ************************************************
//
//Management of Leipzig applicator 
//
#ifndef BrachyDetectorConstructionLeipzig_H
#define BrachyDetectorConstructionLeipzig_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"


class G4LogicalVolume;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Colour;
class G4VPhysicalVolume;
class G4VPhysicalVolume;
class BrachyMaterial;

class BrachyDetectorConstructionLeipzig
{
public:
  BrachyDetectorConstructionLeipzig();
  ~BrachyDetectorConstructionLeipzig(); 
  void  ConstructLeipzig(G4VPhysicalVolume*);

private:
  G4Tubs* capsule ;
  G4LogicalVolume*  capsuleLog;    
  G4VPhysicalVolume* capsulePhys;
 
  G4Sphere* capsuleTip;
  G4LogicalVolume* capsuleTipLog;
  G4VPhysicalVolume* capsuleTipPhys;

  G4Tubs* iridiumCore;
  G4LogicalVolume* iridiumCoreLog;
  G4VPhysicalVolume* iridiumCorePhys;
  
  G4Tubs* applicator1;
  G4LogicalVolume* applicator1Log ;
  G4VPhysicalVolume* applicator1Phys;

  G4Tubs* applicator2;
  G4LogicalVolume* applicator2Log ;
  G4VPhysicalVolume* applicator2Phys;

  BrachyMaterial* pMaterial;   
};

#endif
