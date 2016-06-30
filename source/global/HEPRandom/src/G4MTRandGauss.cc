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

#include <cmath> // for log()

#include "G4MTRandGauss.hh"

// Initialisation of static data
G4ThreadLocal G4bool G4MTRandGauss::set_st = false;
G4ThreadLocal G4double G4MTRandGauss::nextGauss_st = 0.0;

G4MTRandGauss::~G4MTRandGauss()
{
  if ( deleteEngine ) delete localEngine;
}

G4double G4MTRandGauss::operator()()
{
  return fire( defaultMean, defaultStdDev );
}

G4double G4MTRandGauss::operator()( G4double mean, G4double stdDev )
{
  return fire( mean, stdDev );
}

G4double G4MTRandGauss::shoot()
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( getFlag() ) {
    setFlag(false);
    G4double x = getVal();
    return x; 
    // return getVal();
  } 

  G4double r;
  G4double v1,v2,fac,val;
  CLHEP::HepRandomEngine* anEngine = G4MTHepRandom::getTheEngine();

  do {
    v1 = 2.0 * anEngine->flat() - 1.0;
    v2 = 2.0 * anEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = std::sqrt(-2.0*std::log(r)/r);
  val = v1*fac;
  setVal(val);
  setFlag(true);
  return v2*fac;
}

void G4MTRandGauss::shootArray( const G4int size, G4double* vect,
                            G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(mean,stdDev);
}

G4double G4MTRandGauss::shoot( CLHEP::HepRandomEngine* anEngine )
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( getFlag() ) {
    setFlag(false);
    return getVal();
  }

  G4double r;
  G4double v1,v2,fac,val;

  do {
    v1 = 2.0 * anEngine->flat() - 1.0;
    v2 = 2.0 * anEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = std::sqrt( -2.0*std::log(r)/r);
  val = v1*fac;
  setVal(val);
  setFlag(true);
  return v2*fac;
}

void G4MTRandGauss::shootArray( CLHEP::HepRandomEngine* anEngine,
                            const G4int size, G4double* vect,
                            G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = shoot(anEngine,mean,stdDev);
}

G4double G4MTRandGauss::normal()
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( set ) {
    set = false;
    return nextGauss;
  }

  G4double r;
  G4double v1,v2,fac,val;

  do {
    v1 = 2.0 * localEngine->flat() - 1.0;
    v2 = 2.0 * localEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = std::sqrt(-2.0*std::log(r)/r);
  val = v1*fac;
  nextGauss = val;
  set = true;
  return v2*fac;
}

void G4MTRandGauss::fireArray( const G4int size, G4double* vect)
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( defaultMean, defaultStdDev );
}

void G4MTRandGauss::fireArray( const G4int size, G4double* vect,
                           G4double mean, G4double stdDev )
{
   for (G4int i=0; i<size; ++i)
     vect[i] = fire( mean, stdDev );
}

#endif
