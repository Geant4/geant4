// Rich advanced example for Geant4
// RichTbComponent.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbComponent_h
#define RichTbComponent_h 1

#include "globals.hh"
#include "g4std/vector"
#include "RichTbMaterial.hh"
#include "RichTbHall.hh"
#include "RichTbRunConfig.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "AerogelTypeSpec.hh"
#include "G4LogicalBorderSurface.hh"
class RichTbComponent {

public:

   RichTbComponent();
   RichTbComponent(RichTbMaterial*,RichTbHall*, RichTbRunConfig*,G4bool);
  virtual ~RichTbComponent();

   G4LogicalVolume* getEnclosureLogicalVolume()
  {return  RichTbEnclosureLVol;}

   G4VPhysicalVolume* getEnclosurePhysicalVolume()
  {return  RichTbEnclosurePVol;}

   G4LogicalVolume* getMirrorLogicalVolume()
  {return  RichTbMirrorLVol;}

   G4VPhysicalVolume* getMirrorPhysicalVolume()
  {return  RichTbMirrorPVol;}

  G4LogicalBorderSurface* getMirrorSurface()
  {return RichTbMirrorBSurf; }

   G4LogicalVolume* getRadFrameLogicalVolume()
  {return  RichTbRadFrameLVol;}

   G4VPhysicalVolume* getRadFramePhysicalVolume()
  {return  RichTbRadFramePVol;}

   G4LogicalVolume* getRadUpWLogicalVolume()
  {return  RichTbRadUpWLVol;}

   G4VPhysicalVolume* getRadUpWPhysicalVolume()
  {return  RichTbRadUpWPVol;}

   G4LogicalVolume* getRadDnWLogicalVolume()
  {return  RichTbRadDnWLVol;}

   G4VPhysicalVolume* getRadDnWPhysicalVolume()
  {return  RichTbRadDnWPVol;}


  G4int getNumAerogelTiles() {return  NumAerogelTiles; }

   G4LogicalVolume* getAgelLogicalVolume(G4int Agelnum)
  {return  RichTbAgelLVol[Agelnum];}

   G4VPhysicalVolume* getAgelPhysicalVolume(G4int Agelnumber)
  {return  RichTbAgelPVol[Agelnumber];}

   G4LogicalVolume* getAgelWrapTopLogicalVolume(G4int Agelnum)
  {return  RichTbAgelWrapTopLVol[Agelnum];}

   G4VPhysicalVolume* getAgelWrapTopPhysicalVolume(G4int Agelnumber)
  {return  RichTbAgelWrapTopPVol[Agelnumber];}

   G4LogicalVolume* getAgelWrapBotLogicalVolume(G4int Agelnum)
  {return  RichTbAgelWrapBotLVol[Agelnum];}

   G4VPhysicalVolume* getAgelWrapBotPhysicalVolume(G4int Agelnumber)
  {return  RichTbAgelWrapBotPVol[Agelnumber];}

   G4LogicalVolume* getFilterLogicalVolume()
  {return  RichTbFilterLVol;}

   G4VPhysicalVolume* getFilterPhysicalVolume()
  {return  RichTbFilterPVol;}


private:

  G4LogicalVolume* RichTbEnclosureLVol;
  G4VPhysicalVolume* RichTbEnclosurePVol;
  G4LogicalBorderSurface* RichTbEnclosureOuterBSurf;
  G4LogicalBorderSurface* RichTbEnclosureInnerBSurf;  
  G4LogicalVolume* RichTbMirrorLVol;
  G4VPhysicalVolume* RichTbMirrorPVol;
  G4LogicalBorderSurface* RichTbMirrorBSurf;
  G4LogicalVolume* RichTbRadFrameLVol;
  G4VPhysicalVolume* RichTbRadFramePVol;
  G4LogicalVolume* RichTbRadUpWLVol;
  G4VPhysicalVolume* RichTbRadUpWPVol;
  G4LogicalVolume* RichTbRadDnWLVol;
  G4VPhysicalVolume* RichTbRadDnWPVol;
  G4int NumAerogelTiles;
  vector<G4LogicalVolume*> RichTbAgelLVol;
  vector<G4VPhysicalVolume*> RichTbAgelPVol;
  vector<G4LogicalVolume*> RichTbAgelWrapTopLVol;
  vector<G4VPhysicalVolume*> RichTbAgelWrapTopPVol;
  vector<G4LogicalVolume*> RichTbAgelWrapBotLVol;
  vector<G4VPhysicalVolume*> RichTbAgelWrapBotPVol;
  G4LogicalVolume* RichTbFilterLVol;
  G4VPhysicalVolume* RichTbFilterPVol;
  G4bool ConstructTrackingGeometrySwitch;

};

#endif 
