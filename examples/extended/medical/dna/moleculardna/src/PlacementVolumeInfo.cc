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
/// file:PlacementVolumeInfo.cc
/// brief:
#include <utility>
#include "PlacementVolumeInfo.hh"
#include "OctreeNode.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PlacementVolumeInfo::PlacementVolumeInfo(OctreeNode* octree,
                                         OctreeNode* histoneoctree,
                                         std::map<G4int, int64_t> chainbps)
  : fpOctree(octree)
  , fpHistoneOctree(histoneoctree)
  , fBasePairsInChain(std::move(chainbps))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t PlacementVolumeInfo::GetPairsOnChain(G4int idx) const
{
  if(fBasePairsInChain.find(idx) == fBasePairsInChain.end()) {
    return 0;
  } else {
    return fBasePairsInChain.at(idx);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int64_t PlacementVolumeInfo::GetTotalBasePairs() const
{
  int64_t count = 0;
  for(auto it : fBasePairsInChain)
  {
    count += it.second;
  }
  return count;
}