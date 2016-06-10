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
//
// $Id:$

// ------------------------------------------------------------

#include "G4UniformRandPool.hh"
#include "globals.hh"

#include <climits>
#include <stdlib.h>
#include <algorithm>
#include <cstring>

// Not aligned memory
//
void create_pool( G4double*& buffer , G4int ps )
{
  buffer = new G4double[ps];
}

void destroy_pool( G4double*& buffer)
{
  delete[] buffer;
}


#if defined(WIN32)
// No bother with WIN
void create_pool_align( G4double*& buffer , G4int ps)
{
  create_pool(buffer,ps);
}
void destroy_pool_align( G4double*& buffer )
{
  destroy_pool(buffer);
}

#else

// Align memory pools
// Assumption is: static_assert(sizeof(G4double)*CHAR_BIT==64)
//
void create_pool_align( G4double*& buffer , G4int ps)
{
//#if defined(__ICC) || (__INTEL_COMPILER)
//  //INTEL optimized way
//   buffer = (G4doulbe*)mm_allign(ps*sizeof(G4double),sizeof(G4double)*CHAR_BIT);
//#else

  // POSIX standard way
  G4int errcode = posix_memalign( (void**) &buffer ,
                                 sizeof(G4double)*CHAR_BIT,
                                 ps*sizeof(G4double));
  if ( errcode != 0 )
  {
    G4Exception("G4UniformRandPool::create_pool_align()",
                "InvalidCondition", FatalException,
                "Cannot allocate aligned buffer");
    return;
  }
//#endif

  return;
}

void destroy_pool_align( G4double*& buffer )
{
  //#if defined(__ICC) || (__INTEL_COMPILER)
  //mm_free(buffer);
  //#else

  free(buffer);

  //#endif
}
#endif

G4UniformRandPool::G4UniformRandPool()
 : size(G4UNIFORMRANDPOOL_DEFAULT_POOLSIZE), buffer(0), currentIdx(0)
{
  if ( sizeof(G4double)*CHAR_BIT==64 )
  {
    create_pool_align(buffer,size);
  }
  else
  {
    create_pool(buffer,size);
  }
  Fill(size);
}

G4UniformRandPool::G4UniformRandPool( G4int siz )
 : size(siz), buffer(0), currentIdx(0)
{
  if ( sizeof(G4double)*CHAR_BIT==64 )
  {
    create_pool_align(buffer,size);
  }
  else
  {
    create_pool(buffer,size);
  }
  Fill(size);

} 

G4UniformRandPool::~G4UniformRandPool()
{
  if ( sizeof(G4double)*CHAR_BIT==64 )
  {
    destroy_pool_align(buffer);
  }
  else
  {
    destroy_pool(buffer);
  }
}

void G4UniformRandPool::Resize(/*PoolSize_t*/ G4int newSize )
{
  if ( newSize != size )
  {
    destroy_pool(buffer);
    create_pool(buffer,newSize);
    size=newSize;
    currentIdx = 0;
  }
  currentIdx = 0;
}

void G4UniformRandPool::GetMany( G4double* rnds , G4int howmany )
{
  assert(rnds!=0 && howmany>0);

  // if ( howmany <= 0 ) return;
  // We generate at max "size" numbers at once, and
  // We do not want to use recursive calls (expensive).
  // We need to deal with the case  howmany>size
  // So:
  // how many times I need to get "size" numbers?

  const G4int maxcycles = howmany/size;

  // This is the rest
  //
  const G4int peel = howmany%size;
  assert(peel<size);

  // Ok from now I will get random numbers in group of  "size"
  // Note that if howmany<size maxcycles == 0
  //
  G4int cycle = 0;

  // Consider the case howmany>size, then maxcycles>=1
  // and we will request at least "size" rng, so
  // let's start with a fresh buffer of numbers if needed
  //
  if ( maxcycles>0 && currentIdx>0 )
  {
    assert(currentIdx<=size);
    Fill(currentIdx);//<size?currentIdx:size);
  }
  for ( ; cycle < maxcycles ; ++cycle )
  {
    // We can use memcpy of std::copy, it turns out that the two are basically 
    // performance-wise equivalent (expected), since in my tests memcpy is a
    // little bit faster, I use that
    // std::copy(buffer,buffer+size,rnds+(cycle*size));
    //
    memcpy(rnds+(cycle*size),buffer,sizeof(G4double)*size );

    // Get a new set of numbers
    //
    Fill(size);  // Now currentIdx is 0 again
  }

  // If maxcycles>0 last think we did was to call Fill(size)
  // so currentIdx == 0
  // and it is guranteed that peel<size, we have enough fresh random numbers
  // but if maxcycles==0 currentIdx can be whatever, let's make sure we have
  // enough fresh numbers
  //
  if (currentIdx + peel >= size)
  {
    Fill(currentIdx<size?currentIdx:size);
  }
  memcpy(rnds+(cycle*size) , buffer+currentIdx , sizeof(G4double)*peel );
  // std::copy(buffer+currentIdx,buffer+(currentIdx+peel), rnds+(cycle*size));

  // Advance index, we are done
  //
  currentIdx+=peel;
  assert(currentIdx<=size);
}

// Static interfaces implementing CLHEP methods

#include "G4Threading.hh"
#include "G4AutoDelete.hh"

namespace
{
  G4ThreadLocal G4UniformRandPool* rndpool = 0;
}

G4double G4UniformRandPool::flat()
{
  if ( rndpool == 0 )
  {
    rndpool = new G4UniformRandPool;
    G4AutoDelete::Register(rndpool);
  }
  return rndpool->GetOne();
}

void G4UniformRandPool::flatArray( G4int howmany, G4double *rnds)
{
  if ( rndpool == 0 )
  {
    rndpool = new G4UniformRandPool;
    G4AutoDelete::Register(rndpool);
  }
  rndpool->GetMany(rnds,(unsigned int)howmany);
}
