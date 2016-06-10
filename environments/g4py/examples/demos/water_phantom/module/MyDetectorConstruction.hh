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
// $Id: MyDetectorConstruction.hh 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   MyDetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef MY_DETECTOR_CONSTRUCTION_H
#define MY_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================
class G4LogicalVolume;
class G4VSensitiveDetector;

class MyDetectorConstruction : public G4VUserDetectorConstruction {
private:
  G4LogicalVolume* scoreVoxel;

public:
  MyDetectorConstruction();
  ~MyDetectorConstruction(); 

  virtual G4VPhysicalVolume* Construct();
  
  void SetSDtoScoreVoxel(G4VSensitiveDetector* asd);

};

#endif
