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
// $Id: ExN01DetectorConstruction.hh,v 1.3 2006-06-04 21:36:34 kmura Exp $
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

