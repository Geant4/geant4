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
// G4NavigationHistory
//
// Class description:
//
// Responsible for maintenance of the history of the path taken through
// the geometrical hierarchy. Principally a utility class for use by the
// G4Navigator.

// 25.07.96 - P.Kent Initial version. Services derived from
//                   requirements of G4Navigator.
// ----------------------------------------------------------------------
#ifndef G4NAVIGATIONHISTORY_HH
#define G4NAVIGATIONHISTORY_HH

#include <assert.h>
#include <vector>
#include <iostream>

#include "geomdefs.hh"
#include "geomwdefs.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationLevel.hh"
#include "G4NavigationHistoryPool.hh"
#include "G4Allocator.hh"

class G4NavigationHistory
{

 public:  // with description

  friend std::ostream&
  operator << (std::ostream& os, const G4NavigationHistory& h);

  G4NavigationHistory();
    // Constructor: sizes history lists & resets histories.

  ~G4NavigationHistory();
    // Destructor.

  G4NavigationHistory(const G4NavigationHistory& h);
    // Copy constructor.

  inline G4NavigationHistory& operator=(const G4NavigationHistory& h);
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

  inline std::size_t GetDepth() const;
    // Returns current history depth.

  inline std::size_t GetMaxDepth() const;
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

  inline void NewLevel(G4VPhysicalVolume* pNewMother,
                       EVolume vType = kNormal,
                       G4int nReplica = -1);
    // Changes navigation level to that of the new mother.

  inline void BackLevel();
    // Back up one level in history: from mother to grandmother.
    // It does not erase history record of current mother.

  inline void BackLevel(G4int n);
    // Back up specified number of levels in history.

  inline void *operator new(std::size_t);
    // Override "new" for "G4Allocator".
  inline void operator delete(void *aHistory);
    // Override "delete" for "G4Allocator".

 private:

  inline void EnlargeHistory();
    // Enlarge history if required: increase size by kHistoryStride.
    // Note that additional history entries are `dirty' (non zero) apart
    // from the volume history.

 private:

  std::vector<G4NavigationLevel>* fNavHistory;
    // Pointer to the vector of navigation levels.

  std::size_t fStackDepth;
    // Depth of stack: effectively depth in geometrical tree.
};

#include "G4NavigationHistory.icc"

#endif
