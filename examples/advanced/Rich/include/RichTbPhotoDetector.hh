// Rich advanced example for Geant4
// RichTbPhotoDetector.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbPhotoDetector_h
#define RichTbPhotoDetector_h 1

#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbComponent.hh"
#include "RichTbHpd.hh"
#include "RichTbRunConfig.hh"
class RichTbPhotoDetector {
public:

  RichTbPhotoDetector();
  RichTbPhotoDetector(RichTbMaterial*, RichTbComponent*,  RichTbRunConfig*, G4bool);
  virtual ~RichTbPhotoDetector();
  G4int getNumberOfHpds()  {return NumOfHpds; }
  RichTbHpd* getRichHPD(G4int HpdNum){ return richHPD[HpdNum]; }
private:

  G4int NumOfHpds;
  RichTbHpd* richHPD[NumberOfHpds];
  G4bool ConstructTrackingGeometrySwitch;
};

#endif 

