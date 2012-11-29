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
// $Id$
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4PhysicsVectorCache
//
// Class description:
//
// This class includes cache data in use by G4PhysicsVector:
// last input value, last output value, last bin location.

// Author:
// 04.05.2010 Hisaya Kurashige
// ---------------------------------------------------------------

#ifndef G4PhysicsVectorCache_h
#define G4PhysicsVectorCache_h 1

#include "globals.hh"
#include "G4Allocator.hh"

class G4PhysicsVectorCache 
{
  public:  

    G4PhysicsVectorCache();
      // Constructor

   ~G4PhysicsVectorCache();
      // Destructor

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    G4double lastEnergy;        // Cache the last input value
    G4double lastValue;         // Cache the last output value   
    size_t lastBin;             // Cache the last bin location
};

#if defined G4GLOB_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4PhysicsVectorCache> aPVCacheAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4PhysicsVectorCache> aPVCacheAllocator;
#endif

inline void* G4PhysicsVectorCache::operator new(size_t)
{
  void* aCache;
  aCache = (void*)aPVCacheAllocator.MallocSingle();
  return aCache;
}

inline void G4PhysicsVectorCache::operator delete(void* aCache)
{
  aPVCacheAllocator.FreeSingle((G4PhysicsVectorCache*)aCache);
}

#endif
