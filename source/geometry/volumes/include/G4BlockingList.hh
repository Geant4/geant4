// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BlockingList.hh,v 1.2 1999-11-11 15:35:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4BlockingList
//
// A utility class responsible for (efficiently) maintaining a List 
// of blocked volume numbers, with rapid `reset' operations.
//
// Member functions:
// 
// G4BlockingList(G4int maxDefault=511,G4int stride=512)
// 
// 	Create empty blocking List of default size and `stride' resize
//      count.
//
// ~G4BlockingList()
//
// Reset()
//
//	Efficiently `Reset' the blocking List, so that no volumes
//      are blocked
//      [Advance tag no and only fully clear List if tag max reached]
//
// FullyReset()
//
//      Clear the blocking List and reset tag value [slow]
//
// Enlarge(G4int nv)
//
//      Enlarges blocking List is current size < nv, in units of stride
//      Clears the new part of the List
//
// G4int Length()
//     
//      Returns the current length of the List. Note a length of 16
//      means volumes of indices between 0 & 15 inclusive may be blocked
// 
// void BlockVolume(G4int v)
//
// 	Block the volume number v. Requires: 0<=v<Length()
//
// G4bool IsBlocked(G4int v)
//
//      Return true if the volume number v is blocked, else false
//      Requires: 0<=v<Length()
//
// Notes:
//
// Implemented via a ValVector of ints: a tag value is used to set
// the indices of blocked volumes. On reset the current tag value is
// increases, so that the ValVector must only be zeroed when the
// numerical range of the tag is used.
//
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
public:
	G4BlockingList(G4int maxDefault=kBlockingListMaxDefault,
                       G4int stride=kBlockingListStride);
	~G4BlockingList();
	void Reset();
	void FullyReset();
	void Enlarge(const G4int nv);
	G4int Length() const;
	void BlockVolume(const G4int v);
	G4bool IsBlocked(const G4int v) const;

private:
// Current blocked volume tag no.
	G4int fBlockTagNo,fStride;		
// Blocked volumes: Elements with indices
// corresponding to blocked volume set to fBlockTagNo
    	G4RWTValVector<G4int> fBlockingList; 
	
};

#include "G4BlockingList.icc"

#endif
