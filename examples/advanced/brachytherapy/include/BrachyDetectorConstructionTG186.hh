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

class BrachyDetectorConstructionTG186
{
public:
  BrachyDetectorConstructionTG186();
  ~BrachyDetectorConstructionTG186();

  void  ConstructTG186(G4VPhysicalVolume*);
  // Model the TG186 reference source

  void  CleanTG186(); 
  // Destroy the TG186 reference source in the experimental set-up

private:   
  G4Tubs* TG186capsule ;
  G4LogicalVolume*  TG186capsuleLog;    
  G4VPhysicalVolume* TG186capsulePhys;
  G4Sphere* TG186capsuleTip;
  G4LogicalVolume* TG186capsuleTipLog;
  G4VPhysicalVolume* TG186capsuleTipPhys;
  G4Tubs* TG186iridiumCore;
  G4LogicalVolume* TG186iridiumCoreLog;
  G4VPhysicalVolume* TG186iridiumCorePhys;
  G4Tubs* TG186cable;
  G4LogicalVolume* TG186cableLog;
  G4VPhysicalVolume* TG186cablePhys;
    
  BrachyMaterial* pMat;    

  G4VisAttributes* TG186simpleCapsuleVisAtt;
  G4VisAttributes*  TG186simpleCapsuleTipVisAtt;
  G4VisAttributes*  TG186simpleIridiumVisAtt;
  G4VisAttributes*  TG186simpleCableVisAtt;
};
#endif








