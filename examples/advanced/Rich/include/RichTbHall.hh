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
