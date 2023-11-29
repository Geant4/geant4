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
#ifndef MOLECULAR_PLACEMENT_VOLUME_INFO_HH
#define MOLECULAR_PLACEMENT_VOLUME_INFO_HH

#include "globals.hh"

#include <map>

class OctreeNode;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PlacementVolumeInfo
{
 public:
  PlacementVolumeInfo() = default;

  PlacementVolumeInfo(OctreeNode*, OctreeNode*, std::map<G4int, int64_t>);

  virtual ~PlacementVolumeInfo() = default;

  inline auto GetOctree() const { return fpOctree; };

  inline auto GetHistoneOctree() const { return fpHistoneOctree; };

  int64_t GetPairsOnChain(G4int idx) const;

  int64_t GetTotalBasePairs() const;

  inline int64_t GetTotalHistones() const { return fNHistones; }

  inline void SetNHistones(G4int histone) { fNHistones += (int64_t) histone; }

  G4int GetNumberOfChains() const
  {
    return fBasePairsInChain.size();
  };

  void SetOctree(OctreeNode* octree) { fpOctree = octree; };

  void SetHistoneOctree(OctreeNode* octree) { fpHistoneOctree = octree; };

  inline void SetPairsOnChain(G4int chain, int64_t pairs)
  {
    fBasePairsInChain[chain] = pairs;
  };

 private:
  int64_t fNHistones          = 0;
  OctreeNode* fpOctree        = nullptr;
  OctreeNode* fpHistoneOctree = nullptr;
  std::map<G4int, int64_t> fBasePairsInChain;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_PLACEMENT_VOLUME_INFO_HH
