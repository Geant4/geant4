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
//    *    BrachyDetectorConstructionI.hh     *
//    *                                      *
//    ****************************************
// Model of the Iodium source
//
#ifndef BrachyDetectorConstructionI_H
#define BrachyDetectorConstructionI_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class BrachyMaterial;
class G4VisAttributes;

class BrachyDetectorConstructionI
{
public:
  BrachyDetectorConstructionI();
  ~BrachyDetectorConstructionI();
  void  ConstructIodine(G4VPhysicalVolume*);// Construct iodium source
  void  CleanIodine();
private:
  G4Tubs* fDefaultTub;
  G4Tubs* fCapsule;
  G4Sphere* fCapsuleTip;
  G4Tubs* fIodineCore;
  G4LogicalVolume* fDefaultTubLog;
  G4LogicalVolume* fCapsuleLog;
  G4LogicalVolume* fCapsuleTipLog;
  G4LogicalVolume* fIodineCoreLog; 
  G4VPhysicalVolume* fDefaultTubPhys; 
  G4VPhysicalVolume* fCapsulePhys;
  G4VPhysicalVolume* fCapsuleTipPhys1;
  G4VPhysicalVolume* fCapsuleTipPhys2;
  G4VPhysicalVolume* fIodineCorePhys;
  G4VisAttributes* fSimpleIodineVisAtt;
  G4VisAttributes* fSimpleCapsuleVisAtt;
  G4VisAttributes* fSimpleCapsuleTipVisAtt;

  BrachyMaterial* fMaterial;   
};
#endif








