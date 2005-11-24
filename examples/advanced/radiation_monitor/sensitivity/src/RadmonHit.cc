//
// File name:     RadmonHit.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHit.cc,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonHit.hh"
#include "G4Allocator.hh"

void *                                          RadmonHit :: operator new(size_t /* size */)
{
 return reinterpret_cast<void *>(RadmonHitsManager::Instance()->HitsAllocator().MallocSingle());
}
 
 
 
void                                            RadmonHit :: operator delete(void * hit)
{
 RadmonHitsManager::Instance()->HitsAllocator().FreeSingle(reinterpret_cast<RadmonHit *>(hit));
}
 
 
 


G4double                                        RadmonHit :: RetrieveById(G4int id) const
{
 return const_cast<RadmonHit *>(this)->hitData[id];
}
 
 
 
void                                            RadmonHit :: StoreById(G4int id, G4double value)
{
 hitData[id]=value;
}
