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
// $Id: CadDetectorConstruction.hh,v 1.1 2002-06-20 10:00:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------

#ifndef CadDetectorConstruction_h
#define CadDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class CadDetectorMessenger;

class CadDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  CadDetectorConstruction();
  ~CadDetectorConstruction();
  
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
  CadDetectorMessenger* detectorMessenger;
};

#endif

