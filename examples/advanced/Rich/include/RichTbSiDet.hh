// Rich advanced example for Geant4
// RichTbSiDet.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbSiDet_h
#define RichTbSiDet_h 1
#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbGeometryParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbSiPixel.hh"
class RichTbSiDet {
public:
  RichTbSiDet();
  virtual ~RichTbSiDet();
  RichTbSiDet(RichTbMaterial*, G4VPhysicalVolume*, G4bool, G4int, G4int);


  G4LogicalVolume* getHpdSiSectLogicalVolume()
  {return RichTbHpdSiSectLVol;}
  G4VPhysicalVolume* getHpdSiSectPhysicalVolume()
  {return RichTbHpdSiSectPVol;}
  G4LogicalVolume* getHpdSectCoatLogicalVolume()
  {return RichTbHpdSectCoatLVol;}
  G4VPhysicalVolume* getHpdSectCoatPhysicalVolume()
  {return RichTbHpdSectCoatPVol;}
  G4int getNumSiPixels(){return  NumHpdSiPix ;}
  G4int getCurSiSectNum(){return CurSiSectNum;}
  RichTbSiPixel* getHpdSiPx(G4int ipixelnum )
  { return hpdSiPixel[ipixelnum]; }
  G4int getIHpdNumber() 
  { return IHpdNumber; }

private:
  G4LogicalVolume* RichTbHpdSiSectLVol;
  G4VPhysicalVolume* RichTbHpdSiSectPVol;
  G4LogicalVolume* RichTbHpdSectCoatLVol;
  G4VPhysicalVolume* RichTbHpdSectCoatPVol;
  G4int NumHpdSiPix;
  G4int CurSiSectNum;
  RichTbSiPixel* hpdSiPixel[NumberOfPadHpdSiPixels];
  G4bool ConstructTrackGeomSwitch;
  G4int IHpdNumber;

};
#endif 
