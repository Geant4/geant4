// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationHistory.hh,v 1.6 2000-11-01 16:51:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationHistory
//
// Class description:
//
// Responsible for maintenance of the history of the path taken through
// the geometrical hierarchy. Principally a utility class for use by the
// G4Navigator.

// History:
//
// 25.07.96 P.Kent Initial version. Services derived from
//                 requirements of G4Navigator.

#ifndef G4NAVIGATIONHISTORY_HH
#define G4NAVIGATIONHISTORY_HH

#include <assert.h>
#include "geomdefs.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationLevel.hh"

#include "g4rw/tvvector.h"
#include "g4rw/tpvector.h"
#include "g4std/iostream"

class G4NavigationHistory
{

 public:  // with description

  friend G4std::ostream&
  operator << (G4std::ostream &os, const G4NavigationHistory &h);

  G4NavigationHistory();
    // Constructor: sizes history lists & resets histories.

  ~G4NavigationHistory();
    // Destructor.

  G4NavigationHistory(const G4NavigationHistory &h);
    // Copy constructor.

  G4NavigationHistory& operator=(const G4NavigationHistory &h);
    // Assignment operator.

  inline void Reset();
    // Resets history. It now does clear most entries.
    // Level 0 is preserved.

  inline void Clear();
    // Clears entries, zeroing transforms, matrices & negating
    // replica history.

  inline void SetFirstEntry(G4VPhysicalVolume* pVol);
    // Setup initial entry in stack: copies through volume transform & matrix.
    // The volume is assumed to be unrotated.

  inline const G4AffineTransform& GetTopTransform() const; 
    // Returns topmost transform.

  inline const G4AffineTransform* GetPtrTopTransform() const;
    // Returns pointer to topmost transform.

  inline G4int GetTopReplicaNo() const;
    // Returns topmost replica no record.

  inline EVolume GetTopVolumeType() const;
    // Returns topmost volume type.

  inline G4VPhysicalVolume* GetTopVolume() const;
    // Returns topmost physical volume pointer.

  inline G4int GetDepth() const;
    // Returns current history depth.

  inline G4int GetMaxDepth() const;
    // Returns current maximum size of history.
    // Note: MaxDepth of 16 mean history entries [0..15] inclusive.

  inline const G4AffineTransform& GetTransform(G4int n) const;
    // Returns specified transformation.

  inline G4int GetReplicaNo(G4int n) const;
    // Returns specified replica no record.

  inline EVolume GetVolumeType(G4int n) const;
    // Returns specified volume type.

  inline G4VPhysicalVolume* GetVolume(G4int n) const;
    // Returns specified physical volume pointer.

  inline void NewLevel(G4VPhysicalVolume *pNewMother,
		       EVolume vType=kNormal,
		       G4int nReplica=-1);
    // Changes navigation level to that of the new mother.

  inline void BackLevel();
    // Back up one level in history: from mother to grandmother.
    // It does not erase history record of current mother.

  inline void BackLevel(G4int n);
    // Back up specified number of levels in history.

 private:

  inline void EnlargeHistory();
    // Enlarge history if required: increase size by kHistoryStride.
    // Note that additional history entries are `dirty' (non zero) apart
    // from the volume history.

 private:

  G4int fStackDepth;
    // Depth of stack: effectively depth in geometrical tree

  G4RWTValVector<G4NavigationLevel>  fNavHistory;

};

#include "G4NavigationHistory.icc"

#endif
