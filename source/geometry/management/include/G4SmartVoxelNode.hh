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
// G4SmartVoxelNode
//
// Class description:
//
// A node in the smart voxel hierarchy - a 'slice' of space along a given
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

// Author: Paul Kent (CERN), 12.07.1995 - Initial version
//         Gabriele Cosmo (CERN), 18.04.2001 - Migrated to STL vector
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELNODE_HH
#define G4SMARTVOXELNODE_HH

#include <vector>

#include "geomwdefs.hh"
#include "G4Types.hh"

using G4SliceVector = std::vector<G4int>;

/**
 * @brief G4SmartVoxelNode defines a node in the smart voxel hierarchy, i.e.
 * a 'slice' of space along a given axis between given minima and maxima.
 * The node is not aware of its position - this information being
 * available/derivable by the node's owner(s), the voxel headers.
 */

class G4SmartVoxelNode
{
  public:

    /**
     * Constructor for G4SmartVoxelNode. It creates an empty node with slice
     * number 'pSlice'; this number is not stored, but used to provide defaults
     * for the minimum and maximum equivalent node numbers.
     *  @param[in] pSlice Max & min equivalent slice numbers for the header.
     */
    inline G4SmartVoxelNode(G4int pSlice = 0);

    /**
     * Default destructor.
     */
    ~G4SmartVoxelNode() = default;
      // Destructor. No actions.

    /**
     * Equality operator.
     */
    G4bool operator == (const G4SmartVoxelNode& v) const;

    /**
     * Returns the contained volume number 'pVolumeNo'.
     * Note: starts from 0 and no bounds checking is performed.
     */
    inline G4int GetVolume(G4int pVolumeNo) const;

    /**
     * Adds the specified volume number 'pVolumeNo' to the contents.
     */
    inline void Insert(G4int pVolumeNo);

    /**
     * Returns the number of volumes inside the node.
     */
    inline std::size_t GetNoContained() const;

    /**
     * Returns the maximum capacity of the buffer.
     */
    inline std::size_t GetCapacity() const;

    /**
     * Reserves memory in the vector of slices according to the specified
     * quantity, relative to the maximum number of slices.
     */
    inline void Reserve(G4int noSlices);

    /**
     * Shrinks the buffer capacity to the actual size to reduce wasted memory.
     */
    inline void Shrink();

    /**
     * Returns the maximum slice (node/header) number with the same contents
     * and with all intermediate slice also having the same contents.
     */
    inline G4int GetMaxEquivalentSliceNo() const;

    /**
     * Sets the maximum slice number (as above).
     */
    inline void SetMaxEquivalentSliceNo(G4int pMax);

    /**
     * Returns the minimum slice (node/header) number with the same contents
     * and with all intermediate nodes also having the same contents.
     */
    inline G4int GetMinEquivalentSliceNo() const;

    /**
     * Sets the minimum slice number (as above).
     */
    inline void SetMinEquivalentSliceNo(G4int pMin);

  private:

    G4int fminEquivalent;
    G4int fmaxEquivalent;
    G4SliceVector fcontents;
};

#include "G4SmartVoxelNode.icc"

#endif
