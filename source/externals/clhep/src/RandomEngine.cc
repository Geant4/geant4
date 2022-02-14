// -*- C++ -*-
//
// ------------------------------------------------------------------------
//                             HEP Random
//                       --- HepRandomEngine ---
//                      class implementation file
// ------------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// ========================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Moved table of seeds to HepRandom: 19th March 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// =======================================================================

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/EngineFactory.h"

#include <iostream>
#include <vector>

//------------------------- HepRandomEngine ------------------------------

namespace CLHEP {

HepRandomEngine::HepRandomEngine() 
: theSeed (19780503L)
, theSeeds(&theSeed)
{ }

HepRandomEngine::~HepRandomEngine() {}

HepRandomEngine::operator double() {
  return flat();
}

HepRandomEngine::operator float() {
  return float( flat() );
}

HepRandomEngine::operator unsigned int() {
  return (unsigned int)( flat() * exponent_bit_32() );
}

bool 
HepRandomEngine::checkFile (std::istream & file, 
  		  	 const std::string & filename, 
  		  	 const std::string & classname, 
		  	 const std::string & methodname) {
  if (!file) {
    std::cerr << "Failure to find or open file " << filename <<
    " in " << classname << "::" << methodname << "()\n";
    return false;
  }  
  return true;
}			     

std::ostream & HepRandomEngine::put (std::ostream & os) const {
  std::cerr << "HepRandomEngine::put called -- no effect!\n";
  return os;
}
std::istream & HepRandomEngine::get (std::istream & is) {
  std::cerr << "HepRandomEngine::get called -- no effect!\n";
  return is;
}

std::string HepRandomEngine::beginTag ( ) { 
  return "HepRandomEngine-begin"; 
}

std::istream & HepRandomEngine::getState ( std::istream & is ) {
  std::cerr << "HepRandomEngine::getState called -- no effect!\n";
  return is;
}

std::vector<unsigned long> HepRandomEngine::put () const {
  std::cerr << "v=HepRandomEngine::put() called -- no data!\n";
  std::vector<unsigned long> v;
  return v;
}
bool HepRandomEngine::get (const std::vector<unsigned long> & ) {
  std::cerr << "HepRandomEngine::get(v) called -- no effect!\n";
  return false;
}
bool HepRandomEngine::getState (const std::vector<unsigned long> & ) {
  std::cerr << "HepRandomEngine::getState(v) called -- no effect!\n";
  return false;
}

HepRandomEngine* HepRandomEngine::newEngine(std::istream& is) {
  return EngineFactory::newEngine(is);
}

HepRandomEngine* 
HepRandomEngine::newEngine(const std::vector<unsigned long> & v) {
  return EngineFactory::newEngine(v);
}

std::ostream & operator<< (std::ostream & os, const HepRandomEngine & e) {
  return e.put(os);
}

std::istream & operator>> (std::istream & is, HepRandomEngine & e) {
  return e.get(is);
}


}  // namespace CLHEP
