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
// Code developed by
// Silvia Pozzi (1), silvia.pozzi@iss.it
// Barbara Caccia (1), barbara.caccia@iss.it
// Carlo Mancini Terracciano (2), carlo.mancini.terracciano@roma1.infn.it
// (1) Istituto Superiore di Sanita' and INFN Roma, Italy
// (2) Univ. La Sapienza and INFN Roma, Italy

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4RotationMatrix.hh>
#include <G4Material.hh>
#include <G4UnionSolid.hh>
#include <G4PSDoseDeposit3D.hh>
#include <G4MultiFunctionalDetector.hh>

#include "DetectorMessenger.hh"

class G4LogicalVolume;
class DetectorMessenger;
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

	G4VPhysicalVolume* Construct() override;

	void SetJaws(G4double value);
	void SetTargetPosition(G4double value);
	void SetPhantomSide(G4double value);
	void SetVoxelSide(G4double value);
	void SetVoxelDepth(G4double value);

	void UpdateGeometry(G4String, G4double);

	G4int GetNumberSideCells() const;
	G4int GetNumberDepthCells() const;
	G4int GetPhantomDepth() const;

	G4double GetFFilterRadius();
	G4double GetFFilterZ();

	G4double GetJaw1X();
	G4double GetTargetPosition();
	G4double GetAccOriginPosition();
	G4double GetVoxelDepthDim();
    
private:

	void ConstructMaterials();
	G4Material* GetMaterial(G4String materialName);

	void ConstructPhantom();
	void ConstructPhantom_spess();
	void ConstructAccelerator();
	void ConstructTarget();
	void ConstructVacuumWindow();
	void ConstructIonizationChamber();
	void ConstructFlatteningFilter();
	void ConstructPrimaryCollimator();
	void ConstructJawsX();
	void ConstructJawsY();
	void PhantomSegmentation(G4LogicalVolume* phantom);
	G4int CheckPhantomSegmentation(G4int nCells);

	DetectorMessenger* detectorMessenger;
	G4VPhysicalVolume* world_phys;
	G4VPhysicalVolume* accWorld_phys;
	G4VPhysicalVolume* boxJaw1X_phys;
	G4VPhysicalVolume* boxJaw2X_phys;
	G4VPhysicalVolume* boxJaw1Y_phys;
	G4VPhysicalVolume* boxJaw2Y_phys;
	G4VPhysicalVolume* phantom_phys;
	G4LogicalVolume*   LVPhantomSens;

	G4Material* mat_Kapton;
	G4Material* mat_XC10;
	G4Material* mat_WNICU;
	G4Material* mat_Ssteel;
	G4ThreeVector accHalfSize;
	G4ThreeVector jaw1XInitialPos;
	G4ThreeVector jaw2XInitialPos;
	G4ThreeVector jaw1YInitialPos;
	G4ThreeVector jaw2YInitialPos;
	G4double jawAperture;
	G4double fieldSide;
	G4double sourceToSkinDistance;
	G4double phantomSideDim, voxelSideDim, voxelDepthDim;
	G4int nSideCells, nDepthCells;
	G4double tubeFFRadius, tubeFFFirstFaceZ;
	G4MultiFunctionalDetector* phantomDetector3D;
	G4VPrimitiveScorer* phantomScorer3D;
};

#endif
