// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BlockingList.hh,v 1.4 2000-04-25 16:15:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4BlockingList
//
// Class description:
//
// A utility class responsible for (efficiently) maintaining a List 
// of blocked volume numbers, with rapid `reset' operations.
//
// Notes:
//
// Implemented via a ValVector of ints: a tag value is used to set
// the indices of blocked volumes. On reset the current tag value is
// increased, so that the ValVector must only be zeroed when the
// numerical range of the tag is used.

// History:
//
// 24.7.96 P.Kent Separated from G4Navigator

#ifndef G4BLOCKINGLIST_HH
#define G4BLOCKINGLIST_HH

#include "globals.hh"
#include "g4rw/tvvector.h"

const G4int kBlockingListMaxDefault = 500; // Block up to 511 daughters
				           // initially
const G4int kBlockingListStride = 128;
const G4int kBlockTagNoMax = 2147483647;   // 2^31-1 maximum tag no may reach

class G4BlockingList
{
  public:  // with description

    G4BlockingList(G4int maxDefault=kBlockingListMaxDefault,
                   G4int stride=kBlockingListStride);
      // Create empty blocking List of default size and `stride' resize count.

    ~G4BlockingList();
      // Destructor. No operations.

    void Reset();
      // Efficiently `Reset' the blocking List, so that no volumes
      // are blocked [Advance tag number and only fully clear List
      // if tag max reached]

    void FullyReset();
      // Clear the blocking List and reset tag value [slow].

    void Enlarge(const G4int nv);
      // Enlarges blocking List if current size < nv, in units of stride.
      // Clears the new part of the List.

    G4int Length() const;
      // Returns the current length of the List. Note a length of 16
      // means volumes of indices between 0 & 15 inclusive may be blocked.

    void BlockVolume(const G4int v);
      // Block the volume number v.
      // Requires: 0<=v<Length().

    G4bool IsBlocked(const G4int v) const;
      // Return true if the volume number v is blocked, else false.
      // Requires: 0<=v<Length().

  private:

    G4int fBlockTagNo, fStride;		
      // Current blocked volume tag number.

    G4RWTValVector<G4int> fBlockingList; 
      // Blocked volumes: Elements with indices
      // corresponding to blocked volume set to fBlockTagNo.

};

#include "G4BlockingList.icc"

#endif
