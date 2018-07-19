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

#include <assert.h>
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/StaticRandomStates.h"
#include "CLHEP/Utility/memory.h"
#include "CLHEP/Utility/thread_local.h"
#include "CLHEP/Utility/use_atomic.h"

// -----------------------------
// Static members initialisation
// -----------------------------

#include "CLHEP/Random/SeedTable.h"

namespace CLHEP {

  namespace {

    struct defaults {

      defaults()
        : theGenerator( &theDefaultGenerator, do_nothing_deleter() )
        , theEngine   ( &theDefaultEngine, do_nothing_deleter() )
      { }

      defaults(defaults const& other) = delete;
      defaults const& operator=(defaults const&) = delete;

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

    private:

      HepRandom theDefaultGenerator;
      MixMaxRng theDefaultEngine;

    public:

      std::shared_ptr<HepRandom      >  theGenerator;
      std::shared_ptr<HepRandomEngine>  theEngine;
    };  // defaults


#ifdef CLHEP_USE_ATOMIC

    // The ThreadSafeDefaultCache is used only by the function named theDefaults.
    // It is a singly linked list that is intended to hold one object of
    // type "defaults" per thread.

    class ThreadSafeDefaultsCache {
    public:

      ThreadSafeDefaultsCache();

      // The destructor deletes the objects of type "defaults"
      ~ThreadSafeDefaultsCache();

      // Creates new objects and adds them to the linked list in a thread safe manner.
      defaults* createNewDefaults();

      // Note that there are no other functions. No erasing or moving or other accessors.

    private:

      class DefaultsNode {
      public:
        DefaultsNode(DefaultsNode* iNext);
        DefaultsNode const* next() const { return next_; }
        void setNext(DefaultsNode* v) { next_ = v; }
        defaults* addressOfDefaults() { return &defaults_; }
      private:
        DefaultsNode* next_;
        defaults defaults_;
      };

      // points to first node in the linked list
      std::atomic<DefaultsNode*> front_;
    };

    ThreadSafeDefaultsCache::ThreadSafeDefaultsCache() :
      front_(nullptr) {
    }

    defaults* ThreadSafeDefaultsCache::createNewDefaults() {
      DefaultsNode* expected = front_.load();
      DefaultsNode* newNode = new DefaultsNode(expected);
      while (!front_.compare_exchange_strong(expected, newNode)) {
        // another thread changed front_ before us so try again
        newNode->setNext(expected);
      }
      return newNode->addressOfDefaults();
    }

    ThreadSafeDefaultsCache::DefaultsNode::DefaultsNode(DefaultsNode* iNext) :
      next_(iNext),
      defaults_() {
    }

    ThreadSafeDefaultsCache::~ThreadSafeDefaultsCache() {
      DefaultsNode const* node = front_.load();
      while (node) {
        DefaultsNode const* next = node->next();
        delete node;
        node = next;
      }
    }

    defaults &  theDefaults()  {

      // We need to have different engines on different threads because
      // the engines are not thread safe. One cannot generate random numbers
      // using the same engine on different threads simultaneously.
      // Originally we had the defaults object itself as a thread local,
      // but that was failing because on Mac OSX there is not full
      // support for thread locals yet.  Objects containing std::shared_ptr
      // in thread local storage were causing failures. So now we create
      // a container of them that is a function static (not thread local)
      // and the thread local contains only a pointer to an object in the
      // container.
      static ThreadSafeDefaultsCache defaultsForAllThreads;

      // A pointer for each thread to defaults object built for each thread.
      static CLHEP_THREAD_LOCAL defaults* theDefaults = defaultsForAllThreads.createNewDefaults();

      return *theDefaults;
    }
#else

    // This version is used with old compilers not supporting atomics.
    // In that case, the code should not be executed in more than one thread.
    defaults &  theDefaults()  {
      static defaults theDefaults;
      return theDefaults;
    }

#endif

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
