// Rich advanced example for Geant4
// RichTbROGeometry.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbROGeometry.hh"

RichTbROGeometry::RichTbROGeometry(RichTbMaterial* RMaterial,
                                   RichTbRunConfig* RConfig)   
                             : G4VReadOutGeometry() { 
  rMaterial=RMaterial;
  ROTbConfig=RConfig;
 
}


RichTbROGeometry::RichTbROGeometry(G4String aString, 
                           RichTbMaterial* RMaterial, 
                           RichTbRunConfig* RConfig)   
                 : G4VReadOutGeometry(aString ) { 

  rMaterial=RMaterial;
  ROTbConfig=RConfig;
  

}

RichTbROGeometry::~RichTbROGeometry() { }

G4VPhysicalVolume* RichTbROGeometry::Build() {

  //Construct the ReadoutGeometry.
  ROTbHall =new RichTbHall(rMaterial);
  ROTbComponent = new RichTbComponent(rMaterial,ROTbHall,ROTbConfig,false);
  ROPhotDet= new  RichTbPhotoDetector(rMaterial,ROTbComponent,ROTbConfig,false);

  return ROTbHall->getRichTbHallPhysicalVolume();

}
