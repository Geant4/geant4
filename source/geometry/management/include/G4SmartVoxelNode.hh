// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SmartVoxelNode.hh,v 1.1 1999-01-07 16:07:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelNode
//
// A node in the smart voxel hierarchy - a `slice' of space along a given
// axis between given minima and maxima. Note that the node is not aware
// of its position - this information being available/derivable by the
// node's owner(s) (voxelheaders).
//
//
// Member functions:
//
// G4SmartVoxelNode(const G4int pSlice=0)
//   Constructor. Create an empty node with slice number pSlice. THis number
//   is not stored, but used to provide defaults for the minimum and maximum
//   equivalent node numbers
// ~G4SmartVoxelNode()
//   Destructor. No actions.
//
// RWTValOrderedVector<G4int>* GetContents()
//   Return ptr to vector of volume no.s in the node. Use with care.
//   Intended for inspection by navigator at tracking time only.
//
// G4int GetVolume(const G4int pVolumeNo) const
//   Return the pVolumeNo'th contained volume. Note: Starts from 0,
//   no bounds checking.
//
// void Insert(G4int pVolumeNo)
//   Add the specified volume no. to the node's contents
//
// G4int GetNoContained() const
//   Returns the number of volumes contained
//
// G4int GetMaxEquivalenSliceNo() const
//   Return the maximum slice (node/header) no with the same contents, and
//   with all intermediate slice also having the same contents
// void SetMaxEquivalentSliceNo(const G4int pMax)
//   Set the maximum slice no (as above)
// G4int GetMinEquivalentSliceNo() const
//   Return the minimum slice (node/header) no with the same contents, and
//   with all intermediate nodes also having the same contents
// void SetMinEquivalentSliceNo(const G4int pMin)
//   Set the maximum slice no (as above)
//
// Member Data:
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   Min and maximum nodes with same contents. Set by constructor
//   and set methods.
// RWTValOrderedVector<G4int>(1) fcontents
//   Vector of no.s of volumes inside the node
//
// History:
// 12.07.95 P.Kent Initial version

#ifndef G4SMARTVOXELNODE_HH
#define G4SMARTVOXELNODE_HH

#include "globals.hh"
#include "voxeldefs.hh"

#include "G4VPhysicalVolume.hh"
#include <rw/tvordvec.h>


typedef RWTValOrderedVector<G4int> G4SliceVector;

class G4SmartVoxelNode
{
public:

// Constructor. Set min and max equivalent nodes to default.
    G4SmartVoxelNode(const G4int pSlice=0) : fminEquivalent(pSlice),
                                       fmaxEquivalent(pSlice)
    {
    }

// Destructor. No actions necessary
    ~G4SmartVoxelNode()
    {
    }

// Access functions for contents

// Return contained volume no pVolumeNo.
// No bounds checking

    G4int GetVolume(const G4int pVolumeNo) const
    {
	return fcontents(pVolumeNo);
    }

// Add the speicifed volume no to the contents
    void Insert(G4int pVolumeNo)
    {
	fcontents.insert(pVolumeNo);
    }

// Return the no of volumes inside the node
    G4int GetNoContained() const
    {
	return fcontents.entries();
    }

// Access functions for min/max equivalent slices (nodes & headers)
    G4int GetMaxEquivalentSliceNo() const
    {
	return fmaxEquivalent;
    }
    void SetMaxEquivalentSliceNo(const G4int pMax)
    {
	fmaxEquivalent=pMax;
    }
    G4int GetMinEquivalentSliceNo() const
    {
	return fminEquivalent;
    }
    void SetMinEquivalentSliceNo(const G4int pMin)
    {
	fminEquivalent=pMin;
    }

    G4bool operator == (const G4SmartVoxelNode& v) const;
private:
    G4int fminEquivalent;
    G4int fmaxEquivalent;
    G4SliceVector fcontents;
};

#endif


