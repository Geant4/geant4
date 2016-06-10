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
// $Id: G4SmartVoxelNode.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// class G4SmartVoxelNode
//
// Class description:
//
// A node in the smart voxel hierarchy - a `slice' of space along a given
// axis between given minima and maxima. Note that the node is not aware
// of its position - this information being available/derivable by the
// node's owner(s) (voxelheaders).
//
// Member Data:
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   - Min and maximum nodes with same contents. Set by constructor
//     and set methods.
// std::vector<G4int> fcontents
//   - Vector of no.s of volumes inside the node.

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 12.07.95 P.Kent  Initial version
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELNODE_HH
#define G4SMARTVOXELNODE_HH

#include "G4Types.hh"
#include <vector>

typedef std::vector<G4int> G4SliceVector;

class G4SmartVoxelNode
{
  public:  // with description

    G4SmartVoxelNode(G4int pSlice=0) : fminEquivalent(pSlice),
                                       fmaxEquivalent(pSlice) {}
      // Constructor. Create an empty node with slice number pSlice.
      // This number is not stored, but used to provide defaults for the
      // minimum and maximum equivalent node numbers.

    ~G4SmartVoxelNode();
      // Destructor. No actions.

    inline G4int GetVolume(G4int pVolumeNo) const;
      // Return contained volume number pVolumeNo.
      // Note: starts from 0 and no bounds checking performed.

    inline void Insert(G4int pVolumeNo);
      // Add the specified volume number to the contents.

    inline G4int GetNoContained() const;
      // Return the number of volumes inside the node.

    inline G4int GetCapacity() const;
      // Return the maximum capacity of the buffer.

    inline void Reserve(G4int noSlices);
      // Reserve memory in the vector of slices according to the specified
      // quantity, relative to the maximum number of slices.

    inline void Shrink();
      // Shrink buffer capacity to actual size to reduce wasted memory.

    inline G4int GetMaxEquivalentSliceNo() const;
      // Return the maximum slice (node/header) number with the same contents,
      // and with all intermediate slice also having the same contents.
    inline void SetMaxEquivalentSliceNo(G4int pMax);
      // Set the maximum slice number (as above).
    inline G4int GetMinEquivalentSliceNo() const;
      // Return the minimum slice (node/header) number with the same contents,
      // and with all intermediate nodes also having the same contents.
    inline void SetMinEquivalentSliceNo(G4int pMin);
      // Set the minimum slice number (as above).

    G4bool operator == (const G4SmartVoxelNode& v) const;
      // Equality operator.

  private:

    G4int fminEquivalent;
    G4int fmaxEquivalent;
    G4SliceVector fcontents;
};

#include "G4SmartVoxelNode.icc"

#endif
