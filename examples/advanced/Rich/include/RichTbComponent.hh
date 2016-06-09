//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Rich advanced example for Geant4
// RichTbComponent.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbComponent_h
#define RichTbComponent_h 1

#include "globals.hh"
#include <vector>
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
  std::vector<G4LogicalVolume*> RichTbAgelLVol;
  std::vector<G4VPhysicalVolume*> RichTbAgelPVol;
  std::vector<G4LogicalVolume*> RichTbAgelWrapTopLVol;
  std::vector<G4VPhysicalVolume*> RichTbAgelWrapTopPVol;
  std::vector<G4LogicalVolume*> RichTbAgelWrapBotLVol;
  std::vector<G4VPhysicalVolume*> RichTbAgelWrapBotPVol;
  G4LogicalVolume* RichTbFilterLVol;
  G4VPhysicalVolume* RichTbFilterPVol;
  G4bool ConstructTrackingGeometrySwitch;

};

#endif 
