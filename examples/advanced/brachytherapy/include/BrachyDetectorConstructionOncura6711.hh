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
//    Code developed by:
//    S.Guatelli, D. Cutajar, A. Le
// 
//
//    *******************************************
//    *                                          *
//    *  BrachyDetectorConstructionOncura6711.hh *
//    *                                          *
//    ********************************************
//
//    Management of the iodine source
//

#ifndef BrachyDetectorConstruction6711_H
#define BrachyDetectorConstruction6711_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4VPhysicalVolume;
class G4VisAttributes;

class BrachyDetectorConstructionOncura6711
{
public:
  explicit BrachyDetectorConstructionOncura6711();
  ~BrachyDetectorConstructionOncura6711()=default;

  void  ConstructOncura6711(G4VPhysicalVolume*);
  // Model the Oncura 6711 reference source

  void  CleanOncura6711(); 
  // Destroy the Oncura6711 reference source in the experimental set-up

private:   
  G4Tubs* fOncuraCapsule ;
  G4LogicalVolume*  fOncuraCapsuleLog;    
  G4VPhysicalVolume* fOncuraCapsulePhys;
  G4Sphere* fOncuraCapsuleTip1;
  G4LogicalVolume* fOncuraCapsuleTip1Log;
  G4VPhysicalVolume* fOncuraCapsuleTip1Phys;
  G4Sphere* fOncuraCapsuleTip2;
  G4LogicalVolume* fOncuraCapsuleTip2Log;
  G4VPhysicalVolume* fOncuraCapsuleTip2Phys;
  G4Tubs* fOncuraAirGap;
  G4LogicalVolume* fOncuraAirGapLog;
  G4VPhysicalVolume* fOncuraAirGapPhys;
  G4Tubs* fOncuraSilverCore;
  G4LogicalVolume* fOncuraSilverCoreLog;
  G4VPhysicalVolume* fOncuraSilverCorePhys;
  G4VisAttributes* fOncuraCapsuleShellVisAtt;
  G4VisAttributes* fOncuraCapsuleTipVisAtt;
  G4VisAttributes* fOncuraSilverCoreVisAtt;
};
#endif








