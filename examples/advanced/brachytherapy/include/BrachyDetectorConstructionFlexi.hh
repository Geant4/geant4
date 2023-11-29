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
// Code developed by: S. Guatelli, D. Cutajar, J. Poder, 
//Centre For Medical Radiation Physics, University of Wollongong
//
//
//    ******************************************
//    *                                        *
//    *    BrachyDetectorConstructionFlexi.hh  *
//    *                                        *
//    ******************************************
//
// 
//

#ifndef BrachyDetectorConstructionFlexi_H
#define BrachyDetectorConstructionFlexi_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4Tubs;
class G4Cons;
class G4VPhysicalVolume;
class G4VisAttributes;

class BrachyDetectorConstructionFlexi
{
public:
   explicit BrachyDetectorConstructionFlexi();
  ~BrachyDetectorConstructionFlexi()=default;

  void  ConstructFlexi(G4VPhysicalVolume*);
  // Model the Flexi iridium source

  void  CleanFlexi(); 
  // Clean the Iridium source in the experimental set-up

private:   
  G4Tubs* fSteelShell;
  G4LogicalVolume* fLogicalSteelShell;    
  G4VPhysicalVolume* fPhysicalSteelShell;

  G4Tubs* fAirGap;
  G4LogicalVolume* fLogicalAirGap;
  G4VPhysicalVolume* fPhysicalAirGap;

  G4Tubs* fEnd1SteelShell;
  G4LogicalVolume* fLogicalEnd1SteelShell;
  G4VPhysicalVolume* fPhysicalEnd1SteelShell;

  G4Cons* fEnd2SteelShell;
  G4LogicalVolume* fLogicalEnd2SteelShell;
  G4VPhysicalVolume* fPhysicalEnd2SteelShell;

  G4Tubs* fCable;
  G4LogicalVolume* fLogicalCable;
  G4VPhysicalVolume* fPhysicalCable;
    
  G4Tubs* fIridiumCore;
  G4LogicalVolume* fLogicalIridiumCore;
  G4VPhysicalVolume* fPhysicalIridiumCore;

  G4VisAttributes* fSteelAttributes;
  G4VisAttributes* fEndAttributes;
  G4VisAttributes* fSimpleIridiumVisAtt;
};
#endif








