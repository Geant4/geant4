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
// $Id: BrachyDetectorConstructionLeipzig.hh,v 1.2 2002-11-18 15:18:36 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//    ************************************************
//    *                                              *
//    *      BrachyDetectorConstructionLeipzig.hh    *
//    *                                              *
//    ************************************************
//
//Management of Leipzig Applicator
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


 private:

   G4VisAttributes* simpleCapsuleVisAtt;
   G4VisAttributes* simpleCapsuleTipVisAtt;
public:
  
  void  ConstructLeipzig(G4VPhysicalVolume*);

 private:G4Tubs* Capsule ;
  G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
  G4VPhysicalVolume* CapsulePhys;
 
  G4Sphere* CapsuleTip;
  G4LogicalVolume* CapsuleTipLog;
  G4VPhysicalVolume* CapsuleTipPhys;

  G4Tubs* IridiumCore;
  G4LogicalVolume* IridiumCoreLog;
  G4VPhysicalVolume* IridiumCorePhys;
  
  G4Tubs* Appl1;
  G4LogicalVolume* Appl1Log ;
  G4VPhysicalVolume* Appl1Phys;

  G4Tubs* Appl2;
  G4LogicalVolume* Appl2Log ;
  G4VPhysicalVolume* Appl2Phys;

  BrachyMaterial* pMat;   
  };

#endif
