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
// Author: P. Arce
// History: 30.11.07  First version
//*******************************************************
//
// DicomDetectorConstruction.hh :
//	- Start the building of the geometry
//	- Initialisation of materials
//      - Creation of the world 
//	- Reading of the DICOM data
//*******************************************************

#ifndef DicomDetectorConstruction_h
#define DicomDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DicomPhantomZSliceHeader.hh"

class G4Material;
class G4Box;
class G4LogicalVolume;

class DicomDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DicomDetectorConstruction();
  ~DicomDetectorConstruction();

  G4VPhysicalVolume* Construct();
  // trigger the construction of the geometry

protected:
  void InitialisationOfMaterials();
  // create the original materials

  void ReadPhantomData();
  // read the DICOM files describing the phantom

  void ReadPhantomDataFile(const G4String& fname);
  // read one of the DICOM files describing the phantom (usually one per Z slice). Build a DicomPhantomZSliceHeader for each file

  void MergeZSliceHeaders();
  // merge the slice headers of all the files

  G4Material* BuildMaterialWithChangingDensity( const G4Material* origMate, float density, G4String newMateName );
  // build a new material if the density of the voxel is different to the other voxels

  void ConstructPhantomContainer();
  virtual void ConstructPhantom() = 0;
  // construct the phantom volumes. This method should be implemented for each of the derived classes

  void SetScorer(G4LogicalVolume* voxel_logic);

protected:
  G4Material* air;

  // World ...
  G4Box* world_solid;
  G4LogicalVolume* world_logic;
  G4VPhysicalVolume* world_phys;

  G4Box* container_solid;
  G4LogicalVolume* container_logic;
  G4VPhysicalVolume* container_phys;

  G4int fNoFiles; // number of DICOM files
  std::vector<G4Material*> fOriginalMaterials;  // list of original materials 
  std::vector<G4Material*> fMaterials;  // list of new materials created to distinguish different density voxels that have the same original materials
  size_t* fMateIDs; // index of material of each voxel
 //unsigned int* fMateIDs; // index of material of each voxel

  std::map<G4int,G4double> fDensityDiffs; // Density difference to distinguish material for each original material (by index)
 
  std::vector<DicomPhantomZSliceHeader*> fZSliceHeaders; // list of z slice header (one per DICOM files)
  DicomPhantomZSliceHeader* fZSliceHeaderMerged; // z slice header resulted from merging all z slice headers

  G4int nVoxelX, nVoxelY, nVoxelZ;
  G4double voxelHalfDimX,  voxelHalfDimY, voxelHalfDimZ;
};

#endif

