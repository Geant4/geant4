// Rich advanced example for Geant4
// RichTbHpd.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbHpd_h
#define RichTbHpd_h 1
#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbGeometryParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbSiDet.hh"
#include "RichTbSD.hh"
class RichTbHpd {

public:

  RichTbHpd();
  virtual ~RichTbHpd();
  RichTbHpd(RichTbMaterial*, G4VPhysicalVolume * ,G4int, G4int,G4bool);

  G4LogicalVolume* getHpdLogicalVolume()
  {return RichTbHpdLVol;}
  G4VPhysicalVolume* getHpdPhysicalVolume()
  {return RichTbHpdPVol;}
  G4LogicalVolume* getHpdEnvelopeTubeGrandLogicalVolume()
  {return RichTbHpdEnvelopeTubeGrandLVol;}
  G4VPhysicalVolume* getHpdEnvelopeTubeGrandPhysicalVolume()
  {return RichTbHpdEnvelopeTubeGrandPVol;}
  G4LogicalVolume* getHpdQuartzWLogicalVolume()
  {return RichTbHpdQuartzWLVol;}
  G4VPhysicalVolume* getHpdQuartzWPhysicalVolume()
  {return RichTbHpdQuartzWPVol;}
  G4LogicalVolume* getHpdPhCathodeLogicalVolume()
  {return RichTbPhCathodeLVol;}
  G4VPhysicalVolume* getHpdPhCathodePhysicalVolume()
  {return RichTbPhCathodePVol;}
  G4int getHpdNum(){return HpdNumber;}
  G4int getNumSectInHpd() {return  NumSectInHpd;}
  RichTbSiDet* getRichHpdSiDetSect(G4int isect)
  {return richSiliconDetectorSect[isect];}

private:
  G4int HpdNumber;
  RichTbSD* aHpdSD;
  G4LogicalVolume* RichTbHpdLVol;
  G4VPhysicalVolume* RichTbHpdPVol;
  G4LogicalVolume* RichTbHpdEnvelopeTubeGrandLVol;
  G4VPhysicalVolume* RichTbHpdEnvelopeTubeGrandPVol;
  G4LogicalVolume* RichTbHpdQuartzWLVol;
  G4VPhysicalVolume* RichTbHpdQuartzWPVol;
  G4LogicalVolume* RichTbPhCathodeLVol;
  G4VPhysicalVolume* RichTbPhCathodePVol;
  RichTbSiDet* richSiliconDetectorSect[NumberOfSiDetSectors];
  G4int NumSectInHpd;
  G4bool ConstructTrackGeomSwitch;
};
#endif 




