// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NavigationHistory.hh,v 1.5 2000-04-25 16:15:03 gcosmo Exp $
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
#include "globals.hh"

#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationLevel.hh"

#include "g4rw/tvvector.h"
#include "g4rw/tpvector.h"
#include "g4std/iostream"

const G4int kHistoryMax=15;	// Default max size of history
const G4int kHistoryStride=16;   // History increase stride

class G4NavigationHistory
{
 public:  // with description

  friend G4std::ostream& operator << (G4std::ostream &os,const G4NavigationHistory &h);

  G4NavigationHistory();
  G4NavigationHistory(const G4NavigationHistory &h);
  ~G4NavigationHistory();
  
  void Reset();
  void Clear();
  
  void SetFirstEntry(G4VPhysicalVolume* pVol);
  
  const G4AffineTransform& GetTopTransform() const; 

  const G4AffineTransform* GetPtrTopTransform() const;
  G4int GetTopReplicaNo() const;
  EVolume GetTopVolumeType() const;
  G4VPhysicalVolume* GetTopVolume() const;
  
  G4int GetDepth() const;
  G4int GetMaxDepth() const;
  const G4AffineTransform& GetTransform(const G4int n) const;

  G4int GetReplicaNo(const G4int n) const;
  EVolume GetVolumeType(const G4int n) const;
  G4VPhysicalVolume* GetVolume(const G4int n) const;
  
  void NewLevel(G4VPhysicalVolume *pNewMother,
		EVolume vType=kNormal,
		G4int nReplica=-1);
  void BackLevel();
  void BackLevel(G4int n);
  
  G4NavigationHistory& operator=(const G4NavigationHistory &h);

 private:
  void EnlargeHistory();

  G4int fStackDepth;
    // Depth of stack: effectively depth in geometrical tree

  G4RWTValVector<G4NavigationLevel>  fNavHistory;
  // G4RWTPtrVector<G4NavigationLevel>  fNavHistory;

};

#include "G4NavigationHistory.icc"

#endif
