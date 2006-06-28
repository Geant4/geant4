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
// File name:     RadmonHit.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHit.cc,v 1.2 2006-06-28 13:57:09 gunter Exp $
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
