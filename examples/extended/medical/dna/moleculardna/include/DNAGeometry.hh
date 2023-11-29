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
/// file: DNAGeometry.hh
/// brief: This file handls MolecularDNA Geometry

#ifndef MOLECULAR_DNA_GEOMETRY_HH
#define MOLECULAR_DNA_GEOMETRY_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VTouchable.hh"

#include "UtilityFunctions.hh"
#include "MoleculeList.hh"
#include "DamageModel.hh"
#include "ChromosomeMapper.hh"
#include "DNAWorld.hh"
#include "G4VDNAMolecularGeometry.hh"

#include <vector>
#include <map>
#include <unordered_map>
#include <array>

class G4Material;

class G4MolecularConfiguration;

class G4VisAttributes;

class OctreeNode;

class DNAGeometryMessenger;

class PlacementVolumeInfo;

using placement = std::tuple<G4LogicalVolume*, int64_t, int64_t, G4ThreeVector,
                             G4ThreeVector>;  // dousatsu
// G4int, G4int, G4ThreeVector, G4ThreeVector> placement;//ORG

// Holder to describe the global co-ordinates of each DNA chain
// in a given placement volume
// Format << array of global chain indices for local indices 0, 1, 2, 3>
//         < array of start indices for the base pairs along each chain,
//           given by GLOBAL chain index>
//         < array of end indices for the base pairs along each chain,
//           given by GLOBAL chain index>
//         < Boolean to indicate if the strand ID should be flipped due
//           to a 180 degree rotation of the DNA >
//         < Boolean to indicate if volume has actually been placed >
using placement_transform =
  std::tuple<std::array<int, 8>, std::array<int64_t, 8>, std::array<int64_t, 8>,
             G4bool, G4bool>;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct molecule_t
{
  G4String fname;
  G4String fshape;
  int64_t fchain_id;
  int64_t fstrand_id;
  int64_t fbase_idx;
  G4ThreeVector fposition;
  G4ThreeVector frotation;
  G4ThreeVector fsize;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct compareLVByName
{
  bool operator()(const G4LogicalVolume* lhs, const G4LogicalVolume* rhs) const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAGeometry : public G4VDNAMolecularGeometry
{
 public:
  DNAGeometry();

  ~DNAGeometry() override;

  void BuildDNA(G4LogicalVolume*);

  inline DNAWorld* GetDNAWorld();

  void AddVoxelFile(const G4String& vname, const G4String& fname,
                    const G4bool& twist)
  {
    fVoxelNames[vname] = fname;
    fVoxelTwist[vname] = twist;
  }

  void SetFractalFilename(const G4String& fname)
  {
    fFractalCurveFile = fname;
    FillVoxelVectors();
  }

  void SetFractalAnglesAsPi(const G4bool& value) { fAnglesAsPi = value; }

  void EnableCustomMoleculeSizes(const G4bool& value)
  {
    fEnableCustomMoleculeSizes = value;
  }

  void AddChangeMoleculeSize(const G4String& name, const G4ThreeVector& size);

  inline void SetVoxelSideLength(const G4ThreeVector& length)
  { fVoxelSideLength = length; }

  inline void SetFractalScaling(const G4ThreeVector& length)
  { fFractalScaling = length; }

  OctreeNode* GetTopOctreeNode(G4LogicalVolume*) const;

  const PlacementVolumeInfo* GetPVInfo(const G4LogicalVolume*) const;

  inline auto GetChromosomeMapper() const
  { return fpChromosomeMapper; };

  inline void SetOverlaps(G4bool overlap)
  { fCheckOverlaps = overlap; };

  inline G4bool GetOverlaps() const { return fCheckOverlaps; };

  inline void SetDrawCellVolumes(G4bool b) { fDrawCellVolumes = b; };

  inline G4bool GetDrawCellVolumes() const { return fDrawCellVolumes; };

  inline void SetVerbosity(G4int i) { fVerbosity = i; };

  inline G4int GetVerbosity() const { return fVerbosity; };

  inline void SetSmartless(G4int i) { fSmartless = i; };

  inline G4int GetSmartless() const { return fSmartless; };

  inline void SetMediumMaterial(G4Material* mat) { fpMediumMaterial = mat; };

  inline G4LogicalVolume* GetDNAChemVolumePointer() const;

  inline G4LogicalVolume* GetDNAPhysicsVolumePointer() const;

  G4LogicalVolume* GetMatchingChemVolume(G4LogicalVolume*) const;

  G4LogicalVolume* GetMatchingPhysVolume(G4LogicalVolume*) const;

  void PrintOctreeStats();

  inline void SetDirectInteractionRange(G4double r)
  { fDirectInteractionRange = r; };

  void SetRadicalKillDistance(G4double r) { fRadicalKillDistance = r; };

  void SetHistoneScav(G4bool hf) { fUseHistoneScav = hf; };

  inline G4double GetDirectInteractionRange() const
  {
    return fDirectInteractionRange;
  };

  G4double GetRadicalKillDistance() const
  {
    return fRadicalKillDistance;
  };

  const DamageModel* GetDamageModel() { return fpDamageModel; };

  // Tests
  void ChromosomeTest();

  void BasePairIndexTest();

  void UniqueIDTest();

  int64_t GetGlobalUniqueIDTest(  // dousatsu
    int64_t, int64_t, int64_t, int64_t, int64_t) const;   // dousatsu

 private:
  DNAGeometryMessenger* fpMessenger;
  G4String fFractalCurveFile;
  G4bool fAnglesAsPi = false;
  G4bool fEnableCustomMoleculeSizes = false;
  G4bool fCheckOverlaps = false;
  G4bool fDrawCellVolumes = false;
  G4int fVerbosity = 0;
  G4bool fUseHistoneScav = true;
  // vis attrs to draw cells
  G4VisAttributes* fpDrawCellAttrs;
  // Materials
  G4Material* fpMediumMaterial;
  G4Material *fpSugarMaterial, *fpPhosphateMaterial;
  G4Material *fpGuanineMaterial, *fpAdenineMaterial, *fpThymineMaterial;
  G4Material *fpCytosineMaterial, *fpHistoneMaterial;

  // Placement of molecules
  void FillParameterisedSpace();

  G4VPhysicalVolume* PlacePhosphate(G4LogicalVolume*, const molecule_t&,
                                    const molecule_t&);

  G4VPhysicalVolume* PlaceSugar(G4LogicalVolume*, const molecule_t&,
                                const molecule_t&);

  G4VPhysicalVolume* PlaceBase(G4LogicalVolume*, const molecule_t&,
                               const molecule_t&, const molecule_t&,
                               const molecule_t&);

  G4VPhysicalVolume* PlaceHistone(G4LogicalVolume*,  // dousatsu
                                  const molecule_t&);

  // defining voxels
  std::pair<G4LogicalVolume*, PlacementVolumeInfo*> LoadVoxelVolume(
    const G4String&, const G4String&);

  void FillVoxelVectors();

  G4LogicalVolume *fpDNAVolumePhys{}, *fpDNAVolumeChem{};
  std::map<const G4LogicalVolume*, G4LogicalVolume*> fPhysToChem;
  std::map<const G4LogicalVolume*, G4LogicalVolume*> fChemToPhys;
  std::vector<G4ThreeVector> fVoxelPositions;
  std::vector<G4ThreeVector> fVoxelRotations;
  std::vector<G4String> fVoxelTypes;
  std::vector<G4int> fVoxelIndices;
  // std::vector<G4LogicalVolume *> fVoxelLogicalVolumes;
  G4ThreeVector fVoxelSideLength = G4ThreeVector(75 * nm, 75 * nm, 75 * nm);
  G4ThreeVector fFractalScaling = G4ThreeVector(1, 1, 1);
  std::map<G4String, G4String> fVoxelNames;  // For these, the first el.
  std::map<G4String, G4bool> fVoxelTwist;    // is also lv->GetName()

  std::map<const G4LogicalVolume*, PlacementVolumeInfo*> fInfoMap;
  std::map<G4String, G4ThreeVector> fMoleculeSizes;
  G4double fSmartless = 2.0;
  G4double fRadicalKillDistance = 10 * nm,
           fDirectInteractionRange = 6 * angstrom;
  G4double fHistoneInteractionRadius = 25 * angstrom;

  ChromosomeMapper* fpChromosomeMapper;
  DamageModel* fpDamageModel;
  DNAWorld* fpDNAPhysicsWorld;

  // Members for placement transformations
 public:
  inline int64_t GetNumberOfPlacements() const
  {
    return (int64_t) fPlacementTransformations.size();
  }

  G4int GetGlobalChain(G4int vol_idx, G4int local_chain) const
  {
    return std::get<0>(fPlacementTransformations[vol_idx])[local_chain];
  };

  G4int GetLocalChain(G4int vol_idx, G4int global_chain) const;

  inline int64_t GetStartIdx(int64_t vol_idx,
                                    int64_t global_chn) const  // dousatsu
  {
    return std::get<1>(fPlacementTransformations[vol_idx])[global_chn];
  };

  inline int64_t GetEndIdx(int64_t vol_idx,
                                  int64_t global_chn) const  // dousatsu
  {
    return std::get<2>(fPlacementTransformations[vol_idx])[global_chn];
  };

  int64_t GetMaxBPIdx() const;

  inline G4int GetNumberOfChains() const
  {
    return std::get<0>(fPlacementTransformations[0]).size();
  };

  inline G4bool GetStrandsFlipped(G4int vol_idx) const
  {
    return std::get<3>(fPlacementTransformations[vol_idx]);
  };

 private:
  std::vector<placement_transform> fPlacementTransformations;

  void AddNewPlacement(const G4LogicalVolume*, std::array<int, 8>, G4bool,
                       G4bool);

  void AddFourChainPlacement(std::vector<placement>::iterator,
                             std::vector<placement>::iterator, G4bool);

  void AddSingleChainPlacement(std::vector<placement>::iterator,
                               std::vector<placement>::iterator, G4bool);

  // Lookup methods for chemistry handling
 public:
  // Unique ID generators
  inline int64_t GetGlobalPairID(G4int, G4int, int64_t) const;

  int64_t GetGlobalUniqueID(G4VPhysicalVolume*,  // dousatsu
                            const G4VTouchable*) const;

  // And reconstructors from unique ID
  int64_t GetBasePairFromUniqueID(int64_t) const;

  static inline molecule GetMoleculeFromUniqueID(int64_t);

  G4Material* GetMaterialFromUniqueID(int64_t) const;

  G4int GetChainIDFromUniqueID(int64_t) const;

  G4int GetStrandIDFromUniqueID(int64_t) const;

  G4int GetPlacementIndexFromUniqueID(int64_t) const;

  void FindNearbyMolecules(const G4LogicalVolume* motherLogical,
                           const G4ThreeVector& localPosition,
                           std::vector<G4VPhysicalVolume*>& out,
                           double searchRange) override;

  G4bool IsInsideHistone(const G4LogicalVolume* motherLogical,
                         const G4ThreeVector& localPosition) const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline int64_t DNAGeometry::GetGlobalPairID(G4int place_idx, G4int chain_idx,
                                            int64_t base_idx) const
{
  return base_idx + this->GetStartIdx(place_idx, chain_idx);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline molecule DNAGeometry::GetMoleculeFromUniqueID(int64_t idx)
{
  return (molecule)(idx % ::molecule::Count);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4LogicalVolume* DNAGeometry::GetDNAPhysicsVolumePointer() const
{
  return fpDNAVolumePhys;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline DNAWorld* DNAGeometry::GetDNAWorld() { return fpDNAPhysicsWorld; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4LogicalVolume* DNAGeometry::GetDNAChemVolumePointer() const
{
  return fpDNAVolumeChem;
}
#endif  // MOLECULAR_DNA_GEOMETRY_HH
