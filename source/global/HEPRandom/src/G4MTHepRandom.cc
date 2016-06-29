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

#include <CLHEP/Random/StaticRandomStates.h>
#include <CLHEP/Random/JamesRandom.h>

#include "G4MTHepRandom.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

G4ThreadLocal CLHEP::HepRandomEngine* G4MTHepRandom::theEngine = 0;
G4ThreadLocal G4MTHepRandom* G4MTHepRandom::theGenerator = 0;
G4ThreadLocal G4int G4MTHepRandom::isActive  = 0 ;

//---------------------------- HepRandom ---------------------------------

G4MTHepRandom::G4MTHepRandom()
: deleteEngine(false)
{
  createInstance();
}

G4MTHepRandom::G4MTHepRandom(G4long seed)
: deleteEngine(false)
{
  createInstance();
  setTheSeed(seed);
}

G4MTHepRandom::G4MTHepRandom(CLHEP::HepRandomEngine & algorithm)
: deleteEngine(false)
{
  theGenerator = this;
  theEngine = &algorithm;
  isActive = 1;
}

G4MTHepRandom::G4MTHepRandom(CLHEP::HepRandomEngine * algorithm)
: deleteEngine(true)
{
  createInstance();
  theEngine = algorithm;
}

G4MTHepRandom::~G4MTHepRandom()
{
  if ( deleteEngine )  { delete theEngine; theEngine = 0; }
  //if ( theGenerator )  { delete theGenerator; theGenerator = 0; }
}

G4double G4MTHepRandom::flat()
{
  return theEngine->flat();
}

void G4MTHepRandom::flatArray(const G4int size, G4double* vect)
{
  theEngine->flatArray(size,vect);
}

G4double G4MTHepRandom::operator()()
{
  return flat();
}

std::ostream & operator<< (std::ostream & os, const G4MTHepRandom & dist)
{
  return dist.put(os);
}

std::istream & operator>> (std::istream & is, G4MTHepRandom & dist)
{
  return dist.get(is);
}

std::ostream & G4MTHepRandom::put(std::ostream & os) const
{
  return os;
}
std::istream & G4MTHepRandom::get(std::istream & is)
{
  return is;
}

// --------------------------
// Static methods definitions
// --------------------------

void G4MTHepRandom::setTheSeed(G4long seed, G4int lux)
{
  createInstance();
  theEngine->setSeed(seed,lux);
}

G4long G4MTHepRandom::getTheSeed()
{
  createInstance();
  return theEngine->getSeed();
}

void G4MTHepRandom::setTheSeeds(const G4long* seeds, G4int aux)
{
  createInstance();
  theEngine->setSeeds(seeds,aux);
}

const G4long* G4MTHepRandom::getTheSeeds ()
{
  createInstance();
  return theEngine->getSeeds();
}

G4MTHepRandom * G4MTHepRandom::getTheGenerator()
{
  if (!isActive)  { isActive = G4MTHepRandom::createInstanceOnce(); }
  return theGenerator;
}

CLHEP::HepRandomEngine * G4MTHepRandom::getTheEngine()
{
  if (!isActive)  { isActive = G4MTHepRandom::createInstanceOnce(); }
  return theEngine;
}

void G4MTHepRandom::setTheEngine (CLHEP::HepRandomEngine* theNewEngine)
{
  theEngine = theNewEngine;
}

void G4MTHepRandom::saveEngineStatus( const char filename[] )
{
  createInstance();
  theEngine->saveStatus( filename );
}  

void G4MTHepRandom::restoreEngineStatus( const char filename[] )
{
  createInstance();
  theEngine->restoreStatus( filename );
}  

std::ostream& G4MTHepRandom::saveFullState ( std::ostream & os )
{
  os << *getTheEngine();
  return os;
}

std::istream& G4MTHepRandom::restoreFullState ( std::istream & is )
{
  is >> *getTheEngine();
  return is;
}

std::ostream& G4MTHepRandom::saveStaticRandomStates ( std::ostream & os )
{
  return CLHEP::StaticRandomStates::save(os);
}

std::istream& G4MTHepRandom::restoreStaticRandomStates ( std::istream & is )
{
  return CLHEP::StaticRandomStates::restore(is);
}

void G4MTHepRandom::showEngineStatus()
{
  createInstance();
  theEngine->showStatus();
}  

G4int G4MTHepRandom::createInstance()
{
  if (!isActive)  { isActive = G4MTHepRandom::createInstanceOnce(); }
  if (theGenerator) return 1;  // should always be true
  return 0;
}

namespace {
    G4Mutex jamesrandom = G4MUTEX_INITIALIZER;
}

G4int G4MTHepRandom::createInstanceOnce()
{
  if (isActive) return isActive;
  isActive = 1;


  static G4ThreadLocal CLHEP::HepJamesRandom *defaultEngine = 0;
  if (!defaultEngine)  {
      G4AutoLock l(&jamesrandom);
      defaultEngine = new CLHEP::HepJamesRandom;
  }
 
  if ( !theEngine )  { theEngine = defaultEngine; }
  if ( !theGenerator )  { theGenerator = new G4MTHepRandom( theEngine ); }

  return 1;
}

#endif
