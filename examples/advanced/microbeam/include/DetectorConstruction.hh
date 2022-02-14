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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "CellParameterisation.hh"
#include "EMField.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  explicit DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;
  
  void DefineMaterials();
  
  void ConstructSDandField() override;
     
  void SetMassNucleus(G4double mN){ fMassNucleus = mN;}
  G4double GetMassNucleus(){return fMassNucleus;}          

  void SetMassCytoplasm(G4double mC){ fMassCytoplasm = mC;}
  G4double GetMassCytoplasm(){return fMassCytoplasm;}          

  void SetNbOfPixelsInPhantom(G4int nP){ fNbOfPixelsInPhantom = nP;}
  G4int GetNbOfPixelsInPhantom(){return fNbOfPixelsInPhantom;}          

  // Returns the logical volumes

  G4LogicalVolume* GetLogicalCollDetYoke() {return fLogicYoke2;};
  G4LogicalVolume* GetLogicalIsobutane() {return fLogicBoiteIso;};
  G4LogicalVolume* GetLogicalCollDetGap4() {return fLogic4Gap;};
  G4LogicalVolume* GetLogicalPolyprop() {return fLogicBoite3;};
  G4LogicalVolume* GetLogicalKgm() {return fLogicKgm;};
  
  G4Material * GetNucleusMaterial1() {return  fNucleusMaterial1;};
  G4Material * GetNucleusMaterial2() {return  fNucleusMaterial2;};
  G4Material * GetNucleusMaterial3() {return  fNucleusMaterial3;};
  G4Material * GetCytoplasmMaterial1() {return  fCytoplasmMaterial1;};
  G4Material * GetCytoplasmMaterial2() {return  fCytoplasmMaterial2;};
  G4Material * GetCytoplasmMaterial3() {return  fCytoplasmMaterial3;};
  
  const CellParameterisation * GetCellParameterisation() const 
        {return fMyCellParameterisation;};
        
  G4VPhysicalVolume* ConstructLine();     
  
private:
  
  CellParameterisation * fMyCellParameterisation;
  
  // Materials
  G4Material* fDefaultMaterial;
  G4Material* fCollimatorMaterial;
  G4Material* fBoiteMaterial;
  G4Material* fCathodeMaterial;
  G4Material* fVerreMaterial;
  G4Material* fVerre2Material;
  G4Material* fKgmMaterial;
  G4Material* fBoite2Material;
  G4Material* fBoite3Material;
  G4Material* fNucleusMaterial1;
  G4Material* fCytoplasmMaterial1;
  G4Material* fNucleusMaterial2;
  G4Material* fCytoplasmMaterial2;
  G4Material* fNucleusMaterial3;
  G4Material* fCytoplasmMaterial3;
  
// Volumes

  G4VPhysicalVolume* fPhysiWorld;
  G4LogicalVolume*   fLogicWorld;  
  G4Box*             fSolidWorld;
  
  G4VPhysicalVolume* fPhysiVol;
  G4LogicalVolume*   fLogicVol;  
  G4Box*             fSolidVol;
  
  G4VPhysicalVolume* fPhysiBoite;
  G4LogicalVolume*   fLogicBoite;  
  G4Box*             fSolidBoite;
  
  G4VPhysicalVolume* fPhysiYoke1;
  G4LogicalVolume*   fLogicYoke1;  
  G4Box*             fSolidYoke1;

  G4VPhysicalVolume* fPhysi1Gap;
  G4LogicalVolume*   fLogic1Gap;  
  G4Cons*            fSolid1Gap; 

  G4VPhysicalVolume* fPhysi2Gap;
  G4LogicalVolume*   fLogic2Gap;  
  G4Cons*            fSolid2Gap; 
  
  G4VPhysicalVolume* fPhysi3Gap;
  G4LogicalVolume*   fLogic3Gap;  
  G4Cons*            fSolid3Gap; 
  
  G4VPhysicalVolume* fPhysiYoke2;
  G4LogicalVolume*   fLogicYoke2;  
  G4Box*             fSolidYoke2;
  
  G4VPhysicalVolume* fPhysi4Gap;
  G4LogicalVolume*   fLogic4Gap;  
  G4Cons*            fSolid4Gap; 

  G4VPhysicalVolume* fPhysi5Gap;
  G4LogicalVolume*   fLogic5Gap;  
  G4Cons*            fSolid5Gap; 
    
  G4VPhysicalVolume* fPhysiBoiteIso;
  G4LogicalVolume*   fLogicBoiteIso;  
  G4Box*             fSolidBoiteIso;
  
  G4VPhysicalVolume* fPhysiCathode;
  G4LogicalVolume*   fLogicCathode;  
  G4Box*             fSolidCathode;
  
  G4VPhysicalVolume* fPhysiIso;
  G4LogicalVolume*   fLogicIso;  
  G4Box*             fSolidIso;
  
  G4VPhysicalVolume* fPhysiVerre;
  G4LogicalVolume*   fLogicVerre;  
  G4Box*             fSolidVerre;

  G4VPhysicalVolume* fPhysiBoite2;
  G4LogicalVolume*   fLogicBoite2;  
  G4Box*             fSolidBoite2;

  G4VPhysicalVolume* fPhysiBoite3;
  G4LogicalVolume*   fLogicBoite3;  
  G4Box*             fSolidBoite3;
  
  G4VPhysicalVolume* fPhysiKgm;
  G4LogicalVolume*   fLogicKgm;  
  G4Box*             fSolidKgm;

  G4VPhysicalVolume* fPhysiVerre2;
  G4LogicalVolume*   fLogicVerre2;  
  G4Box*             fSolidVerre2;

  // CELL PHANTOM

  G4VPhysicalVolume* fPhysiPhantom;
  G4LogicalVolume*   fLogicPhantom;  
  G4Box*             fSolidPhantom; 

  G4double fMassNucleus;
  G4double fMassCytoplasm;

  G4double fDensityPhantom;
  G4double fDensityNucleus;
  G4double fDensityCytoplasm;
    
  G4double fWorldSizeXY;
  G4double fWorldSizeZ;
  G4double fCollObjSizeXY;
  G4double fCollObjSizeZ;

  G4double fCiblePositionX;
  G4double fCiblePositionY;
  G4double fCiblePositionZ;
  
  G4double fLineAngle;
  
  G4int    fNbOfPixelsInPhantom;
};

#endif
