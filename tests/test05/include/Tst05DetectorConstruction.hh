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
// $Id: Tst05DetectorConstruction.hh,v 1.4 2001-07-11 10:09:38 gunter Exp $
// ------------------------------------------------------------

#ifndef Tst05DetectorConstruction_h
#define Tst05DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class Tst05DetectorMessenger;

class Tst05DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst05DetectorConstruction();
  ~Tst05DetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  void SelectDetector(G4String val);
  
private:
  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;
  
  G4double calBox_x;
  G4double calBox_y;
  G4double calBox_z;
  G4double rotAngle;
  G4double calPos;

  G4double trackerRadius;
  G4double trackerHight;
  G4double trackerPos;

  G4int detectorChoice;
  Tst05DetectorMessenger* detectorMessenger;
};

#endif

