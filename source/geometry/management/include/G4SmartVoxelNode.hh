//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SmartVoxelNode.hh,v 1.7 2002-04-19 08:20:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// G4std::vector<G4int> fcontents
//   - Vector of no.s of volumes inside the node.

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 12.07.95 P.Kent  Initial version

#ifndef G4SMARTVOXELNODE_HH
#define G4SMARTVOXELNODE_HH

#include "globals.hh"
#include "voxeldefs.hh"

#include "G4VPhysicalVolume.hh"
#include "g4std/vector"

typedef G4std::vector<G4int> G4SliceVector;

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

    G4int GetVolume(G4int pVolumeNo) const;
      // Return contained volume number pVolumeNo.
      // Note: starts from 0 and no bounds checking performed.

    void Insert(G4int pVolumeNo);
      // Add the specified volume number to the contents.

    G4int GetNoContained() const;
      // Return the number of volumes inside the node.

    G4int GetMaxEquivalentSliceNo() const;
      // Return the maximum slice (node/header) number with the same contents,
      // and with all intermediate slice also having the same contents.
    void SetMaxEquivalentSliceNo(G4int pMax);
      // Set the maximum slice number (as above).
    G4int GetMinEquivalentSliceNo() const;
      // Return the minimum slice (node/header) number with the same contents,
      // and with all intermediate nodes also having the same contents.
    void SetMinEquivalentSliceNo(G4int pMin);
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
