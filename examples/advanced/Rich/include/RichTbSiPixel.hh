// Rich advanced example for Geant4
// RichTbSiPixel.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbSiPixel_h
#define RichTbSiPixel_h 1
#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbGeometryParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
class RichTbSiPixel {
public:
  RichTbSiPixel();
  virtual ~RichTbSiPixel();
  RichTbSiPixel(RichTbMaterial*, G4VPhysicalVolume*, G4int, G4int, G4int);

  G4LogicalVolume* getHpdSiPixelLogicalVolume()
  {return RichTbHpdSiPixelLVol;}
  G4VPhysicalVolume* getHpdSiPixelPhysicalVolume()
  {return RichTbHpdSiPixelPVol;}
  G4int getICurSectorNumber() {return ICurSectorNumber; }
  G4int getHpdSiPixelnum() {return HpdSiPixelnum; }
  G4int getICurHpdNumber() {return  ICurHpdNumber;}
private:
  G4LogicalVolume* RichTbHpdSiPixelLVol;
  G4VPhysicalVolume* RichTbHpdSiPixelPVol;
  G4int ICurSectorNumber;
  G4int HpdSiPixelnum;
  G4int ICurHpdNumber;
  G4bool PixelIsAlive;
};
#endif 
