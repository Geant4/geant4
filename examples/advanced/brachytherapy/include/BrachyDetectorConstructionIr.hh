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
// $Id: BrachyDetectorConstructionIr.hh,v 1.4 2006-05-12 13:49:29 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionIr.hh   *
//    *                                      *
//    ****************************************
//
// Management of the Iridium source
//

#ifndef BrachyDetectorConstructionIr_H
#define BrachyDetectorConstructionIr_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4VPhysicalVolume;
class BrachyMaterial;
class G4VisAttributes;

class BrachyDetectorConstructionIr 
{
public:
  BrachyDetectorConstructionIr();
  ~BrachyDetectorConstructionIr();

  void  ConstructIridium(G4VPhysicalVolume*);
  // Model the Iridum source

  void  CleanIridium(); 
  // Destroy the Iridium source in the experimental set-up

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
    
  BrachyMaterial* pMat;    

  G4VisAttributes* simpleCapsuleVisAtt;
  G4VisAttributes*  simpleCapsuleTipVisAtt;
  G4VisAttributes*  simpleIridiumVisAtt;
};
#endif








