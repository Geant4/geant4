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
// Class G4SmartVoxelNode implementation
//
// 18.04.01, G.Cosmo - Migrated to STL vector
// 12.07.95, P.Kent - Initial version
// --------------------------------------------------------------------

#include "G4SmartVoxelNode.hh"

// --------------------------------------------------------------------
G4Allocator<G4SmartVoxelNode>*& aNodeAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4SmartVoxelNode>* _instance = nullptr;
  return _instance;
}

// --------------------------------------------------------------------
// Return true if contents equal
//
// Preconditions:
//
// Node contents were entered in the same order
//
G4bool G4SmartVoxelNode::operator == (const G4SmartVoxelNode& v) const
{
  std::size_t maxNode = GetNoContained();
  if (maxNode == v.GetNoContained())
  {
    for (std::size_t node=0; node<maxNode; ++node)
    {
      if (GetVolume((G4int)node) != v.GetVolume((G4int)node))
      {
        return false;
      }
    }
    return true;
  }
  return false;
}
