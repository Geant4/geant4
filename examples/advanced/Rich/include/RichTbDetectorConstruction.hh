// Rich advanced example for Geant4
// RichTbDetectorConstruction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////

#ifndef RichTbDetectorConstruction_h
#define RichTbDetectorConstruction_h 1

#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbHall.hh"
#include "RichTbComponent.hh"
#include "RichTbPhotoDetector.hh"
#include "RichTbGraphics.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbRunConfig.hh"
#include "RichTbROGeometry.hh"
#include "RichTbAnalysisManager.hh"
class RichTbDetectorConstruction:
      public G4VUserDetectorConstruction{

  public: 
  RichTbDetectorConstruction();

  RichTbDetectorConstruction(RichTbRunConfig*);
  virtual ~RichTbDetectorConstruction();
  G4VPhysicalVolume* Construct();
  RichTbMaterial* getRichTbMaterial()
  {return  rMaterial; }
  RichTbHall* getRichTbHall() 
  {return rTbHall; }
  RichTbComponent* getRichTbComponent()
  {return  rTbComponent; }
  RichTbPhotoDetector* getRichTbPhotoDetector()
  {return rTbPhotoDetector; }
  RichTbGraphics* getRichTbGraphcis()
  {return  rTbGraphics; }
  RichTbROGeometry* getROGeometry() 
  {return  rTbROGeom; }
  private:

  RichTbRunConfig* runConfiguration;
  RichTbMaterial* rMaterial;
  RichTbHall* rTbHall;
  RichTbComponent* rTbComponent;
  RichTbPhotoDetector* rTbPhotoDetector;
  RichTbGraphics* rTbGraphics;
  RichTbROGeometry* rTbROGeom;


};

#endif 




