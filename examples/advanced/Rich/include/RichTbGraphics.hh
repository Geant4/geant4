// Rich advanced example for Geant4
// RichTbGraphics.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbGraphics_h
#define RichTbGraphics_h 1

#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbHall.hh"
#include "RichTbComponent.hh"
#include "RichTbPhotoDetector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "RichTbROGeometry.hh"

class RichTbGraphics {

public:

   RichTbGraphics();
  RichTbGraphics(RichTbHall*, RichTbComponent*, RichTbPhotoDetector*,
                 RichTbROGeometry*, RichTbRunConfig* );
  virtual ~RichTbGraphics();

private:

};

#endif 
