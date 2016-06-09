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
