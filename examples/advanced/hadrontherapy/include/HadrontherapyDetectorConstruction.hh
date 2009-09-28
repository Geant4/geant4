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
// $Id: HadrontherapyDetectorConstruction.hh; 
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyDetectorConstruction_H
#define HadrontherapyDetectorConstruction_H 1

#include "globals.hh"
#include "G4VisAttributes.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class HadrontherapyDetectorROGeometry;
class HadrontherapyDetectorMessenger;
class HadrontherapyDetectorSD;

class HadrontherapyDetectorConstruction 
{
public:

  HadrontherapyDetectorConstruction(G4VPhysicalVolume*);

  ~HadrontherapyDetectorConstruction();


private: 

  void ConstructPhantom();
  void ConstructDetector();
  void ConstructSensitiveDetector();
 
  //  G4VisAttributes* redWire;
  
public: 

    G4LogicalVolume* GetDetectorLogicalVolume(){ return detectorLogicalVolume;}

  // This method allows to set phantom geometry
  // void SetPhantomDimensionAndPosition(G4double );
  // This method allows to set detector geometry
  // void SetDetectorDimensionAndPosition(G4double );

    G4double ComputeVoxelSize() {return detectorSizeX/numberOfVoxelsAlongX;};
  // Returns the size of the voxel along the X axis
 
private:

  HadrontherapyDetectorMessenger* detectorMessenger; 

  G4VisAttributes* skyBlue;
  G4VisAttributes* red;

  G4VPhysicalVolume* mother;

  HadrontherapyDetectorSD* detectorSD; // Pointer to sensitive detector
  HadrontherapyDetectorROGeometry* detectorROGeometry; // Pointer to ROGeometry 


  G4double phantomSizeX; 
  G4double phantomSizeY; 
  G4double phantomSizeZ;

  G4double detectorSizeX; 
  G4double detectorSizeY; 
  G4double detectorSizeZ;


  G4VPhysicalVolume* phantomPhysicalVolume;

  G4LogicalVolume*   detectorLogicalVolume;
  G4VPhysicalVolume* detectorPhysicalVolume;


  G4int numberOfVoxelsAlongX; 
  G4int numberOfVoxelsAlongY;
  G4int numberOfVoxelsAlongZ;  
};
#endif
