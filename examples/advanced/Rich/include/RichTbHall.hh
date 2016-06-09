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
// Rich advanced example for Geant4
// RichTbHall.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbHall_h
#define RichTbHall_h 1

#include "globals.hh"
#include "RichTbMaterial.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
class RichTbHall {

public:
  RichTbHall();
  RichTbHall(RichTbMaterial*);
  virtual  ~RichTbHall();
  G4LogicalVolume* getRichTbHallLogicalVolume()
  {return RichTbHallLVol;}
  G4VPhysicalVolume* getRichTbHallPhysicalVolume()
  {return RichTbHallPVol;}
  
private:
  
  G4LogicalVolume* RichTbHallLVol;
  G4VPhysicalVolume* RichTbHallPVol;


};

#endif 
