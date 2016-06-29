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
#if __clang__
  #if ((defined(G4MULTITHREADED) && !defined(G4USE_STD11)) || \
      !__has_feature(cxx_thread_local)) || !__has_feature(c_atomic)
    #define CLANG_NOSTDTLS
  #endif
#endif

#if (defined(G4MULTITHREADED) && \
    (!defined(G4USE_STD11) || (defined(CLANG_NOSTDTLS) || defined(__INTEL_COMPILER))))

#include "G4MTRandFlat.hh"

const G4int G4MTRandFlat::MSBBits= 15;
const unsigned long G4MTRandFlat::MSB= 1ul<<G4MTRandFlat::MSBBits;
G4ThreadLocal unsigned long G4MTRandFlat::staticRandomInt= 0;
G4ThreadLocal unsigned long G4MTRandFlat::staticFirstUnusedBit= 0;

G4MTRandFlat::~G4MTRandFlat()
{
  if ( deleteEngine ) delete localEngine;
}

G4double G4MTRandFlat::operator()()
{
  return fire( defaultA, defaultB );
}

G4double G4MTRandFlat::operator()( G4double w )
{
  return fire( w );
}

G4double G4MTRandFlat::operator()( G4double a, G4double b )
{
  return fire( a, b );
}

G4double G4MTRandFlat::shoot()
{
  return G4MTHepRandom::getTheEngine()->flat();
}

void G4MTRandFlat::shootArray(const G4int size, G4double* vect)
{
  G4MTHepRandom::getTheEngine()->flatArray(size,vect);
}

void G4MTRandFlat::shootArray( const G4int size, G4double* vect,
                               G4double lx, G4double dx  )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(lx,dx);
}

void G4MTRandFlat::shootArray( CLHEP::HepRandomEngine* anEngine,
                               const G4int size, G4double* vect,
                               G4double lx, G4double dx  )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(anEngine,lx,dx);
}

void G4MTRandFlat::fireArray( const G4int size, G4double* vect)
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( defaultA, defaultB );
}

void G4MTRandFlat::fireArray( const G4int size, G4double* vect,
                              G4double lx, G4double dx  )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( lx, dx );
}

#endif
