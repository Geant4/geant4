// Rich advanced example for Geant4
// RichTbROGeometry.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbROGeometry_h
#define RichTbROGeometry_h 1
#include "G4VReadOutGeometry.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbMaterial.hh"
#include "RichTbHall.hh"
#include "RichTbComponent.hh"
#include "RichTbRunConfig.hh"
#include "RichTbPhotoDetector.hh"

class RichTbROGeometry : public G4VReadOutGeometry {
public:
  RichTbROGeometry(RichTbMaterial*,RichTbRunConfig*);
  RichTbROGeometry(G4String, RichTbMaterial*,RichTbRunConfig*);
  ~RichTbROGeometry();

  G4VPhysicalVolume* Build();

  RichTbHall* getRichTbHall() 
  { return ROTbHall; }
  RichTbComponent* getRichTbComponent() 
  { return ROTbComponent; }
  RichTbPhotoDetector* getROPhotDet()
  { return  ROPhotDet ;}  
private:
  
  RichTbMaterial* rMaterial;
  RichTbRunConfig* ROTbConfig;
  RichTbHall* ROTbHall;
  RichTbComponent* ROTbComponent;
  RichTbPhotoDetector*  ROPhotDet;

};
#endif
