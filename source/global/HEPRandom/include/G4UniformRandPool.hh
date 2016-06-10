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
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Class description:
//
// 

// ------------------------------------------------------------
#ifndef G4UNIFORMRANDPOOL_HH
#define G4UNIFORMRANDPOOL_HH

#include <algorithm>
#include <assert.h>

#include "G4Types.hh"
#include "Randomize.hh"

#define G4UNIFORMRANDPOOL_DEFAULT_POOLSIZE 1024
#define G4UNIFORMRANDPOOL_TINY_POOLSIZE 128
#define G4UNIFORMRANDPOOL_SMALL_POOLSIZE 256
#define G4UNIFORMRANDPOOL_MEDIUM_POOLSIZE 512
#define G4UNIFORMRANDPOOL_LARGE_POOLSIZE 2048
#define G4UNIFORMRANDPOOL_HUGE_POOLSIZE 8192

class G4UniformRandPool
{
  public:

    G4UniformRandPool();
    G4UniformRandPool( /*PoolSize_t&*/ G4int ps );
   ~G4UniformRandPool();

    void Resize( /*PoolSize_t*/ G4int newSize );
    void GetMany( G4double* rnds , G4int howMany );
    inline G4double GetOne();
    inline /*PoolSize_t*/ G4int GetPoolSize() const;

    // These two static methods are used to
    // simulate the calls of CLHEP::HepRandom
    //
    static G4double flat();
    static void flatArray( G4int howmany, G4double* rnds );

  private:

    inline void Fill( G4int howmany );

  private:

    /*PoolSize_t*/ G4int size;
    G4double* buffer;
    G4int currentIdx;
};

inline G4double G4UniformRandPool::GetOne()
{
  // No more available numbers, re-fill
  //
  if ( currentIdx >= /*(unsigned int)*/size )
  {
    Fill(/*(unsigned int)*/size);
  }

  return buffer[currentIdx++];
}

inline /*PoolSize_t*/ G4int G4UniformRandPool::GetPoolSize() const
{
  return size;
}

inline void G4UniformRandPool::Fill( G4int howmany )
{
  assert(howmany>0 && howmany <= size);

  // Fill buffer with random numbers
  //
  G4Random::getTheEngine()->flatArray(howmany,buffer);
  currentIdx = 0;
}

#endif
