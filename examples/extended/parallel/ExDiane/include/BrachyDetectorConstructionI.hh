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
// $Id: BrachyDetectorConstructionI.hh,v 1.2 2004/05/25 08:36:17 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionI.hh     *
//    *                                      *
//    ****************************************
//Author: Susanna Guatelli
//
#ifndef BrachyDetectorConstructionI_H
#define BrachyDetectorConstructionI_H 1

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

class BrachyDetectorConstructionI
{
public:
  BrachyDetectorConstructionI();
  ~BrachyDetectorConstructionI();
  void  ConstructIodium(G4VPhysicalVolume*);// Construct iodium source

private:
  //air tub ...
  G4Tubs* defaultTub; 
  G4LogicalVolume* defaultTubLog;
  G4VPhysicalVolume* defaultTubPhys;
  
  //source titanium capsule ...
  G4Tubs* capsule ;
  G4LogicalVolume*  capsuleLog;   
  G4VPhysicalVolume* capsulePhys;

  G4Sphere* capsuleTip;
  G4LogicalVolume* capsuleTipLog;
  G4VPhysicalVolume* capsuleTipPhys1;
  G4VPhysicalVolume* capsuleTipPhys2;

  //radioactive core ...
  G4Tubs* iodiumCore;
  G4LogicalVolume* iodiumCoreLog;
  G4VPhysicalVolume* iodiumCorePhys;

  // golden marker ...
  G4Tubs* marker;
  G4LogicalVolume* markerLog;
  G4VPhysicalVolume* markerPhys;

  BrachyMaterial* pMaterial;   
};
#endif








