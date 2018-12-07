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
/// \file medical/DICOM/include/DicomPartialDetectorConstruction.hh
/// \brief Definition of the DicomPartialDetectorConstruction class
//

#ifndef DicomPartialDetectorConstruction_h
#define DicomPartialDetectorConstruction_h 1

#include <map>

#include "globals.hh"
#include "DicomDetectorConstruction.hh"

class G4PartialPhantomParameterisation;
class G4LogicalVolume;
class G4Material;

/// Construct a DICOM Geant4 geometry produced from the intersetion
/// of a DICOM file and a volume

class DicomPartialDetectorConstruction : public DicomDetectorConstruction
{

public:

  DicomPartialDetectorConstruction();
  ~DicomPartialDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

private:

  virtual void ReadPhantomData();
  void ReadPhantomDataFile(const G4String& fname);
  void ConstructPhantomContainer();
  virtual void ConstructPhantom();

  void ReadVoxelDensitiesPartial( std::ifstream& fin );

// void ReadVoxelDensitiesPartial( std::ifstream& fin,
//                      std::map< G4int, std::map< G4int, G4int > > ifxmin,
//                      std::map< G4int, std::map< G4int, G4int > > ifxmax );
//
// std::pair<G4double,G4double> ReadVoxelDim(G4int nVoxel, std::ifstream& fin);
//
// void SetScorer(G4LogicalVolume* voxel_logic);
//
// G4Material* BuildMaterialChangingDensity(const G4Material* origMate, 
//                                          G4float density, G4String mateName);

private:

  G4PartialPhantomParameterisation* fPartialPhantomParam;

  std::multimap<G4int,G4int> fFilledIDs;
  std::map< G4int, std::map< G4int, G4int > > fFilledMins;
  std::map< G4int, std::map< G4int, G4int > > fFilledMaxs;
  G4int fNVoxels;
  G4double fDimX, fDimY, fDimZ;
  G4double fOffsetX, fOffsetY, fOffsetZ;
  std::vector<G4Material*> fPhantomMaterials;
};

#endif
