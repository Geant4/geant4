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
class BrachyMaterial;
class G4VisAttributes;

class BrachyDetectorConstructionOncura6711
{
public:
  BrachyDetectorConstructionOncura6711();
  ~BrachyDetectorConstructionOncura6711();

  void  ConstructOncura6711(G4VPhysicalVolume*);
  // Model the Oncura 6711 reference source

  void  CleanOncura6711(); 
  // Destroy the Oncura6711 reference source in the experimental set-up

private:   
  G4Tubs* OncuraCapsule ;
  G4LogicalVolume*  OncuraCapsuleLog;    
  G4VPhysicalVolume* OncuraCapsulePhys;
  G4Sphere* OncuraCapsuleTip1;
  G4LogicalVolume* OncuraCapsuleTip1Log;
  G4VPhysicalVolume* OncuraCapsuleTip1Phys;
  G4Sphere* OncuraCapsuleTip2;
  G4LogicalVolume* OncuraCapsuleTip2Log;
  G4VPhysicalVolume* OncuraCapsuleTip2Phys;
  G4Tubs* OncuraAirGap;
  G4LogicalVolume* OncuraAirGapLog;
  G4VPhysicalVolume* OncuraAirGapPhys;
  G4Tubs* OncuraSilverCore;
  G4LogicalVolume* OncuraSilverCoreLog;
  G4VPhysicalVolume* OncuraSilverCorePhys;

  BrachyMaterial* pMat;    

  G4VisAttributes* OncuraCapsuleShellVisAtt;
  G4VisAttributes*  OncuraCapsuleTipVisAtt;
  G4VisAttributes*  OncuraSilverCoreVisAtt;
};
#endif








