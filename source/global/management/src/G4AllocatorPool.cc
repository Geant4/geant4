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
// $Id: G4AllocatorPool.cc,v 1.3 2005/03/15 19:11:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
// G4AllocatorPool
//
// Implementation file
//
// Author: G.Cosmo, November 2000
//

#include "G4AllocatorPool.hh"

// ************************************************************
// G4AllocatorPool constructor
// ************************************************************
//
G4AllocatorPool::G4AllocatorPool( unsigned int sz )
  : esize(sz<sizeof(G4PoolLink) ? sizeof(G4PoolLink) : sz),
    csize(sz<1024/2-16 ? 1024-16 : sz*10-16),
    chunks(0), head(0), nchunks(0)
{
}

// ************************************************************
// G4AllocatorPool copy constructor
// ************************************************************
//
G4AllocatorPool::G4AllocatorPool(const G4AllocatorPool& right)
  : esize(right.esize), csize(right.csize)
{
  *this = right;
}

// ************************************************************
// G4AllocatorPool operator=
// ************************************************************
//
G4AllocatorPool&
G4AllocatorPool::operator= (const G4AllocatorPool& right)
{
  if (&right == this) { return *this; }
  chunks  = right.chunks;
  head    = right.head;
  nchunks = right.nchunks;
  return *this;
}

// ************************************************************
// G4AllocatorPool destructor
// ************************************************************
//
G4AllocatorPool::~G4AllocatorPool()
{
  Reset();
}

// ************************************************************
// Reset
// ************************************************************
//
void G4AllocatorPool::Reset()
{
  // Free all chunks
  //
  G4PoolChunk* n = chunks;
  G4PoolChunk* p = 0;
  while (n)
  {
    p = n;
    n = n->next;
    delete p;
  }
  head = 0;
  chunks = 0;
  nchunks = 0;
}

// ************************************************************
// Grow
// ************************************************************
//
void G4AllocatorPool::Grow()
{
  // Allocate new chunk, organize it as a linked list of
  // elements of size 'esize'
  //
  G4PoolChunk* n = new G4PoolChunk(csize);
  n->next = chunks;
  chunks = n;
  nchunks++;

  const int nelem = csize/esize;
  char* start = n->mem;
  char* last = &start[(nelem-1)*esize];
  for (char* p=start; p<last; p+=esize)
  {
    reinterpret_cast<G4PoolLink*>(p)->next
      = reinterpret_cast<G4PoolLink*>(p+esize);
  }
  reinterpret_cast<G4PoolLink*>(last)->next = 0;
  head = reinterpret_cast<G4PoolLink*>(start);
}
