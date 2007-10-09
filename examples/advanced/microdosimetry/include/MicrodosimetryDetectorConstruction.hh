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
// -------------------------------------------------------------------
// $Id: MicrodosimetryDetectorConstruction.hh,v 1.2 2007-10-09 08:08:42 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrodosimetryDetectorConstruction_h
#define MicrodosimetryDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"

#include "MicrodosimetryPhantomConfiguration.hh"
#include "MicrodosimetryCellParameterisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrodosimetryDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  MicrodosimetryDetectorConstruction();
  ~MicrodosimetryDetectorConstruction();

  G4VPhysicalVolume* Construct();
     
  void SetMassNucleus(G4float mN){ massNucleus = mN;}
  G4float GetMassNucleus(){return massNucleus;}          

  void SetMassCytoplasm(G4float mC){ massCytoplasm = mC;}
  G4float GetMassCytoplasm(){return massCytoplasm;}          

  void SetNbOfPixelsInPhantom(G4int nP){ nbOfPixelsInPhantom = nP;}
  G4int GetNbOfPixelsInPhantom(){return nbOfPixelsInPhantom;}          

private:

  G4float massPhantom;
  G4float massNucleus;
  G4float massCytoplasm;
  G4double densityPhantom;
  G4double densityNucleus;
  G4double densityCytoplasm;
  G4int nbOfPixelsInPhantom;
    
  G4double      WorldSizeXY;
  G4double      WorldSizeZ;
  G4double      CollObjSizeXY;
  G4double      CollObjSizeZ;

 
// Materials

  G4Material*   defaultMaterial;
  G4Material*   KgmMaterial;

// Volumes

  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  
   G4VPhysicalVolume* physiKgm;
  G4LogicalVolume*   logicKgm;  
  G4Box*             solidKgm;
// CELL

  G4VPhysicalVolume* physiPhantom;
  G4LogicalVolume*   logicPhantom;  
  G4Box*             solidPhantom; 

  MicrodosimetryPhantomConfiguration myMicrodosimetryPhantomConfiguration;
  MicrodosimetryCellParameterisation* phantomParam ; 

  void DefineMaterials();
  G4VPhysicalVolume* ConstructMicrodosimetryLine();     

};

#endif
