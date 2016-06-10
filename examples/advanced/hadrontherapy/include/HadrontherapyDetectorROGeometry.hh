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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

//The detectior is devided in voxels. 
//
//
#ifndef HadrontherapyDetectorROGeometry_h
#define HadrontherapyDetectorROGeometry_h 

//#include "G4VReadOutGeometry.hh"
#include "G4VUserParallelWorld.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

class HadrontherapyDetectorROGeometry : public G4VUserParallelWorld
{
public:
  HadrontherapyDetectorROGeometry(G4String);
  ~HadrontherapyDetectorROGeometry();

  void Initialize(G4ThreeVector detectorPos,
		  G4double detectorDimX,
		  G4double detectorDimY,
		  G4double detectorDimZ,
		  G4int numberOfVoxelsX,
		  G4int numberOfVoxelsY,
		  G4int numberOfVoxelsZ);

  void UpdateROGeometry();

  virtual void Construct();
  virtual void ConstructSD();

private:  
  //Parameters used for the construction
  G4ThreeVector detectorToWorldPosition; 
  G4double detectorSizeX;
  G4double detectorSizeY; 
  G4double detectorSizeZ;

  G4int numberOfVoxelsAlongX;
  G4int numberOfVoxelsAlongY; 
  G4int numberOfVoxelsAlongZ; 
  
  //Solids that are updated on-the fly
  G4Box* RODetector;
  G4Box* RODetectorXDivision;
  G4Box* RODetectorYDivision;
  G4Box* RODetectorZDivision;

  //Logical volumes used for the re-build on-the-fly
  G4LogicalVolume* worldLogical;
  G4LogicalVolume* RODetectorLog;
  G4LogicalVolume* RODetectorXDivisionLog;
  G4LogicalVolume* RODetectorYDivisionLog;
  G4LogicalVolume* RODetectorZDivisionLog;
  G4LogicalVolume* sensitiveLogicalVolume;


  G4bool isBuilt;
  G4bool isInitialized;
};
#endif
