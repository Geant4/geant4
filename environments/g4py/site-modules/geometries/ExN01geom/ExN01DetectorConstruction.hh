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
// $Id: ExN01DetectorConstruction.hh,v 1.4 2006-06-29 15:29:12 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   ExN01DetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EXN01_DETECTOR_CONSTRUCTION_H
#define EXN01_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class ExN01DetectorConstruction : public G4VUserDetectorConstruction {
public:
  
  ExN01DetectorConstruction();
  ~ExN01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  
private:
  // Logical volumes
  //
  G4LogicalVolume* experimentalHall_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* calorimeterBlock_log;
  G4LogicalVolume* calorimeterLayer_log;
  
  // Physical volumes
  //
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* calorimeterLayer_phys;
  G4VPhysicalVolume* calorimeterBlock_phys;
  G4VPhysicalVolume* tracker_phys;

};

#endif

