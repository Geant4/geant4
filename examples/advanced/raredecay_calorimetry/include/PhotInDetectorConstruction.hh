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
//
// $Id: PhotInDetectorConstruction.hh,v 1.3 2006/06/29 16:24:41 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// 

#ifndef PhotInDetectorConstruction_h
#define PhotInDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "PhotInLayerParameterisation.hh"
#include "PhotInGapParameterisation.hh"
#include "PhotInCalorimeterSD.hh"
#include "PhotInConstants.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class PhotInDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  // Constructor ^ Distructor
  PhotInDetectorConstruction(G4double x=1., G4double y=1., G4double z=1.);
  virtual ~PhotInDetectorConstruction();

  // Interface for the Geometry Construction (can be overloaded)
  virtual G4VPhysicalVolume* Construct();

  // Member functions
  void CreateMaterial(G4String materialChoice);
  void UpdateGeometry();
  void PrintCalorParameters() const;
  void SetAbsorberMaterial(G4String materialChoice);     
  void SetGapMaterial(G4String materialChoice);     
  void SetSerialGeometry(G4bool ser);
  void SetNumberOfLayers(G4int ln);
  void SetNumberOfSlabs(G4int sn);
  void SetSamplingFraction(G4double sf);
  //void SetSectionHalfDimensions(G4double xhd,G4double xhd, G4double xhd);
		// Immediate definition of Get functions
  G4double GetHalfZThickness()   { return zHD; }
  G4double GetHalfYWidth()       { return yHD; }
  G4double GetHalfXHight()       { return xHD; }
  G4String GetAbsorberMaterial() const { return absorberMaterial->GetName();}
  G4String GetGapMaterial() const      { return gapMaterial->GetName();};
  G4int GetNumberOfLayers() const      { return numberOfLayers; }
  G4int GetNumberOfSlabs() const       { return numberOfSlabs; }
  G4double GetSamplingFraction() const { return samplingFraction; }
  G4bool IsSerial() const              { return serial; }
     
private:
  void DefineMaterials();

		// ---- Body ----
  G4int                      numberOfLayers;   // #of layers in one section of Calorimeter
  G4int                      numberOfSlabs;    // #of slabs of the active gap in a layer
  G4double                   samplingFraction; // gap_thickness/total_thickness
  G4double                   xHD;              // Dimension of the section along the slab
  G4double                   yHD;              // Dimension of the section perpend the slab
  G4double                   zHD;              // Dimension of the section perpen the layer

  G4bool                     serial;           // A flag of serial or parallel layout

  G4Material*                worldMaterial;    // Material of the WORLD volume (?shape)
  G4Material*                absorberMaterial; // Material of the Absorber (variable)
  G4Material*                gapMaterial;      // Material of the Active Gap (variable)

  G4Box*                     layerSolid;       // Layer of absorber including a  gap
  G4Box*                     slabSolid;        // Active Gap (SD=Sensitive Detector)
  G4LogicalVolume*           calorLogical[PhotInNumSections]; // logVol's for CalorSections
  G4LogicalVolume*           layerLogical[PhotInNumSections]; // logVol's for LayersOfSects
  G4LogicalVolume*           slabLogical[PhotInNumSections];  // logVol's for SlabsOfLayers
  G4VPhysicalVolume*         calorPhysical[PhotInNumSections];//physVol's for CalSections
  G4PVParameterised*         layerPhysical[PhotInNumSections];//physVol's for LayersOfSects
  G4PVParameterised*         slabPhysical[PhotInNumSections]; //physVol's for SlabsOfLayers

  PhotInLayerParameterisation* layerParam;     // Parameterization of layers
  PhotInGapParameterisation*   gapParam;       // Parameterization of gaps
      
};

#endif

