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
// $Id: BrachyDetectorConstructionI.hh,v 1.2 2002-11-18 15:18:35 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionI.hh     *
//    *                                      *
//    ****************************************
//Management of the Iodium source

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


private:
G4VisAttributes*   simpleMarkerVisAtt;
G4VisAttributes* simpleIodiumVisAtt;

G4VisAttributes* simpleCapsuleVisAtt;
G4VisAttributes* simpleCapsuleTipVisAtt;

public:


void  ConstructIodium(G4VPhysicalVolume*);
private:
G4Tubs* DefaultTub;
G4LogicalVolume* DefaultTub_log;
G4VPhysicalVolume*DefaultTub_Phys;

G4Tubs* Capsule ;
G4LogicalVolume*  CapsuleLog;    //pointer to the logical World
G4VPhysicalVolume* CapsulePhys;

G4Sphere* CapsuleTip;
G4LogicalVolume* CapsuleTipLog;
G4VPhysicalVolume* CapsuleTipPhys1;
G4VPhysicalVolume* CapsuleTipPhys2;

G4Tubs* IodiumCore;
G4LogicalVolume* ICoreLog;
G4VPhysicalVolume* ICorePhys;

G4Tubs* Marker;
G4LogicalVolume* MarkerLog;
G4VPhysicalVolume* MarkerPhys;

BrachyMaterial* pMat;   
};



#endif








