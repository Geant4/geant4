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
// $Id$
// -------------------------------------------------------------------

#ifndef MicrobeamDetectorConstruction_h
#define MicrobeamDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
//#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"

#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"

#include "MicrobeamPhantomConfiguration.hh"
#include "MicrobeamCellParameterisation.hh"
#include "MicrobeamEMField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  MicrobeamDetectorConstruction();
  ~MicrobeamDetectorConstruction();

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

  G4double 	CiblePositionX;
  G4double 	CiblePositionY;
  G4double 	CiblePositionZ;
  
  G4double 	lineAngle;
 
// Materials

  G4Material*   defaultMaterial;
  G4Material*   collimatorMaterial;
  G4Material*   BoiteMaterial;
  G4Material*   CathodeMaterial;
  G4Material*   VerreMaterial;
  G4Material*   Verre2Material;
  G4Material*   KgmMaterial;
  G4Material*   Boite2Material;
  G4Material*   Boite3Material;
  G4Material*   nucleusMaterial1;
  G4Material*   cytoplasmMaterial1;
  G4Material*   nucleusMaterial2;
  G4Material*   cytoplasmMaterial2;
  G4Material*   nucleusMaterial3;
  G4Material*   cytoplasmMaterial3;

// Volumes

  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logicWorld;  
  G4Box*             solidWorld;
  
  G4VPhysicalVolume* physiVol;
  G4LogicalVolume*   logicVol;  
  G4Box*             solidVol;
  
  G4VPhysicalVolume* physiBoite;
  G4LogicalVolume*   logicBoite;  
  G4Box*             solidBoite;
  
  G4VPhysicalVolume* physiYoke1;
  G4LogicalVolume*   logicYoke1;  
  G4Box*             solidYoke1;

  G4VPhysicalVolume* physi1Gap;
  G4LogicalVolume*   logic1Gap;  
  G4Cons*            solid1Gap; 

  G4VPhysicalVolume* physi2Gap;
  G4LogicalVolume*   logic2Gap;  
  G4Cons*            solid2Gap; 
  
  G4VPhysicalVolume* physi3Gap;
  G4LogicalVolume*   logic3Gap;  
  G4Cons*            solid3Gap; 
  
  G4VPhysicalVolume* physiYoke2;
  G4LogicalVolume*   logicYoke2;  
  G4Box*             solidYoke2;
  
  G4VPhysicalVolume* physi4Gap;
  G4LogicalVolume*   logic4Gap;  
  G4Cons*            solid4Gap; 

  G4VPhysicalVolume* physi5Gap;
  G4LogicalVolume*   logic5Gap;  
  G4Cons*            solid5Gap; 
    
  G4VPhysicalVolume* physiBoiteIso;
  G4LogicalVolume*   logicBoiteIso;  
  G4Box*             solidBoiteIso;
  
  G4VPhysicalVolume* physiCathode;
  G4LogicalVolume*   logicCathode;  
  G4Box*             solidCathode;
  
  G4VPhysicalVolume* physiIso;
  G4LogicalVolume*   logicIso;  
  G4Box*             solidIso;
  
  G4VPhysicalVolume* physiVerre;
  G4LogicalVolume*   logicVerre;  
  G4Box*             solidVerre;

  G4VPhysicalVolume* physiBoite2;
  G4LogicalVolume*   logicBoite2;  
  G4Box*             solidBoite2;

  G4VPhysicalVolume* physiBoite3;
  G4LogicalVolume*   logicBoite3;  
  G4Box*             solidBoite3;
  
  G4VPhysicalVolume* physiKgm;
  G4LogicalVolume*   logicKgm;  
  G4Box*             solidKgm;

  G4VPhysicalVolume* physiVerre2;
  G4LogicalVolume*   logicVerre2;  
  G4Box*             solidVerre2;
    
// CELL

  G4VPhysicalVolume* physiPhantom;
  G4LogicalVolume*   logicPhantom;  
  G4Box*             solidPhantom; 

  MicrobeamPhantomConfiguration myMicrobeamPhantomConfiguration;
  MicrobeamCellParameterisation* phantomParam ; 
  MicrobeamEMField * Field;

  void DefineMaterials();
  G4VPhysicalVolume* ConstructMicrobeamLine();     

  G4FieldManager *pFieldMgr;
  G4MagIntegratorStepper * pStepper;
  G4EqMagElectricField * pEquation;
  G4MagInt_Driver * pIntgrDriver;
  G4ChordFinder *pChordFinder ;
  G4PropagatorInField *propInField;

};

#endif
