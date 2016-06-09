// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- HepRandom ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - HepRandom defined as singleton, constructors are
//                  kept public for backward compatibility: 27th Feb 1998
//                - Relocated Poisson and Gauss data and simplified
//                  initialisation of static generator: 5th Jan 1999
// =======================================================================

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/StaticRandomStates.h"
#include "CLHEP/Utility/memory.h"

// -----------------------------
// Static members initialisation
// -----------------------------

#include "CLHEP/Random/SeedTable.h"

namespace CLHEP {


namespace {

struct defaults {
  defaults( HepRandom & g, HepJamesRandom & e )
    : theGenerator( &g, do_nothing_deleter() )
    , theEngine   ( &e, do_nothing_deleter() )
  { }

  void  resetEngine( HepRandomEngine * newEngine ) {
    theEngine.reset( newEngine );
  }

  void  resetEngine( HepRandomEngine & newEngine ) {
    theEngine.reset( &newEngine, do_nothing_deleter() );
  }

  bool  ensureInitialized()  {
    assert( theGenerator.get() != 0  && theEngine.get() != 0 );
    return true;
  }

  ~defaults()
  { }

  shared_ptr<HepRandom      >  theGenerator;
  shared_ptr<HepRandomEngine>  theEngine;
};  // defaults

  inline
  defaults &  theDefaults()  {
    static  HepRandom       theDefaultGenerator;
    static  HepJamesRandom  theDefaultEngine;
    static  defaults theDefaults(theDefaultGenerator, theDefaultEngine);
    return theDefaults;
  }

}  // namespace

//---------------------------- HepRandom ---------------------------------

HepRandom::HepRandom()
{ }

HepRandom::HepRandom(long seed)
{
  setTheSeed(seed);
}

HepRandom::HepRandom(HepRandomEngine & algorithm)
{
  theDefaults().resetEngine( algorithm );
}

HepRandom::HepRandom(HepRandomEngine * algorithm)
{
  theDefaults().resetEngine( algorithm );
}

HepRandom::~HepRandom() 
{ }

double HepRandom::flat()
{
  return theDefaults().theEngine->flat();
}

void HepRandom::flatArray(const int size, double* vect)
{
  theDefaults().theEngine->flatArray(size,vect);
}

double HepRandom::operator()() {
  return flat();
}

std::string HepRandom::name() const {return "HepRandom";}
HepRandomEngine & HepRandom::engine() {
  std::cerr << "HepRandom::engine() called -- there is no assigned engine!\n";
  return *theDefaults().theEngine.get();
} 

std::ostream & operator<< (std::ostream & os, const HepRandom & dist) {
  return dist.put(os);
}

std::istream & operator>> (std::istream & is, HepRandom & dist) {
  return dist.get(is);
}

std::ostream & HepRandom::put(std::ostream & os) const {return os;}
std::istream & HepRandom::get(std::istream & is) {return is;}

// --------------------------
// Static methods definitions
// --------------------------

void HepRandom::setTheSeed(long seed, int lux)
{
  theDefaults().theEngine->setSeed(seed,lux);
}

long HepRandom::getTheSeed()
{
  return theDefaults().theEngine->getSeed();
}

void HepRandom::setTheSeeds(const long* seeds, int aux)
{
  theDefaults().theEngine->setSeeds(seeds,aux);
}

const long* HepRandom::getTheSeeds ()
{
  return theDefaults().theEngine->getSeeds();
}

void HepRandom::getTheTableSeeds(long* seeds, int index)
{
  if ((index >= 0) && (index < 215)) {
    seeds[0] = seedTable[index][0];
    seeds[1] = seedTable[index][1];
  }
  else seeds = NULL;
}

HepRandom * HepRandom::getTheGenerator()
{
  return theDefaults().theGenerator.get();
}

HepRandomEngine * HepRandom::getTheEngine()
{
  return theDefaults().theEngine.get();
}

void HepRandom::setTheEngine (HepRandomEngine* theNewEngine)
{
  theDefaults().theEngine.reset( theNewEngine, do_nothing_deleter() );
}

void HepRandom::saveEngineStatus( const char filename[] )
{
  theDefaults().theEngine->saveStatus( filename );
}  

void HepRandom::restoreEngineStatus( const char filename[] )
{
  theDefaults().theEngine->restoreStatus( filename );
}  

std::ostream& HepRandom::saveFullState ( std::ostream & os ) {
  os << *getTheEngine();
  return os;
}

std::istream& HepRandom::restoreFullState ( std::istream & is ) {
  is >> *getTheEngine();
  return is;
}

std::ostream& HepRandom::saveStaticRandomStates ( std::ostream & os ) {
  return StaticRandomStates::save(os);
}

std::istream& HepRandom::restoreStaticRandomStates ( std::istream & is ) {
  return StaticRandomStates::restore(is);
}

void HepRandom::showEngineStatus()
{
  theDefaults().theEngine->showStatus();
}  

int HepRandom::createInstance()
{
  return static_cast<int>( theDefaults().ensureInitialized() );
}

}  // namespace CLHEP
