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
