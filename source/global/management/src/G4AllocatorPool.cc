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
// G4AllocatorPool class implementation
//
// Author: G.Cosmo, November 2000
// --------------------------------------------------------------------

#include "G4AllocatorPool.hh"

// ************************************************************
// G4AllocatorPool constructor
// ************************************************************
//
G4AllocatorPool::G4AllocatorPool(unsigned int sz)
  : esize(sz < sizeof(G4PoolLink) ? sizeof(G4PoolLink) : sz)
  , csize(sz < 1024 / 2 - 16 ? 1024 - 16 : sz * 10 - 16)
{}

// ************************************************************
// G4AllocatorPool operator=
// ************************************************************
//
G4AllocatorPool& G4AllocatorPool::operator=(const G4AllocatorPool& right)
{
  if(&right == this)
  {
    return *this;
  }
  chunks  = right.chunks;
  head    = right.head;
  nchunks = right.nchunks;
  return *this;
}

// ************************************************************
// G4AllocatorPool destructor
// ************************************************************
//
G4AllocatorPool::~G4AllocatorPool() { Reset(); }

// ************************************************************
// Reset
// ************************************************************
//
void G4AllocatorPool::Reset()
{
  // Free all chunks
  //
  G4PoolChunk* n = chunks;
  G4PoolChunk* p = nullptr;
  while(n != nullptr)
  {
    p = n;
    n = n->next;
    delete p;
  }
  head    = nullptr;
  chunks  = nullptr;
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
  auto* n        = new G4PoolChunk(csize);
  n->next        = chunks;
  chunks         = n;
  ++nchunks;

  const int nelem = csize / esize;
  char* start     = n->mem;
  char* last      = &start[(nelem - 1) * esize];
  for(char* p = start; p < last; p += esize)
  {
    reinterpret_cast<G4PoolLink*>(p)->next =
      reinterpret_cast<G4PoolLink*>(p + esize);
  }
  reinterpret_cast<G4PoolLink*>(last)->next = nullptr;
  head = reinterpret_cast<G4PoolLink*>(start);
}
