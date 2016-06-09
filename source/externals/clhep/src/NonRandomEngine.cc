// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                           Hep Random
//                        --- NonRandomEngine ---
//                   class implementation file
// -----------------------------------------------------------------------
// M. Fischler    - Created 9/30/99
//
// M. Fischler    - Modifications to capture sequence as a vector, which
//		    are needed to retain sanity when put and get are involved.
// Mark Fischler  - Methods for distrib. instance save/restore 12/8/04    
// M. Fischler	  - Initialization of all state data (even those parts unused)
//                - at ctor time, to thwart a VC++ i/o bug.
// M. Fischler    - put/get for vectors of ulongs		3/15/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//
//=========================================================================

#include "CLHEP/Random/NonRandomEngine.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Random/DoubConv.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>

namespace CLHEP {

std::string NonRandomEngine::name() const {return "NonRandomEngine";}

NonRandomEngine::NonRandomEngine() : nextHasBeenSet(false), 
				     sequenceHasBeenSet(false),
				     intervalHasBeenSet(false) ,
                                     nextRandom(0.05),
				     nInSeq(0),
				     randomInterval(0.1) { }

NonRandomEngine::~NonRandomEngine() { }


void NonRandomEngine::setNextRandom(double r) {
  nextRandom = r;
  nextHasBeenSet=true;
  return;
}

void NonRandomEngine::setRandomSequence(double* s, int n) {
  sequence.clear();
  for (int i=0; i<n; i++) sequence.push_back(*s++);
  assert (sequence.size() == (unsigned int)n);
  nInSeq = 0;
  sequenceHasBeenSet=true;
  nextHasBeenSet=false;
  return;
}

void NonRandomEngine::setRandomInterval(double x) {
  randomInterval = x;
  intervalHasBeenSet=true;
  return;
}

double NonRandomEngine::flat() {

  if (sequenceHasBeenSet) {
    double v = sequence[nInSeq++];
    if (nInSeq >= sequence.size() ) sequenceHasBeenSet = false;
    return v;
  }

  if ( !nextHasBeenSet ) {
    std::cout 
	<< "Attempt to use NonRandomEngine without setting next random!\n";
    exit(1);
  }

  double a = nextRandom;
  nextHasBeenSet = false;

  if (intervalHasBeenSet) {
    nextRandom += randomInterval;
    if ( nextRandom >= 1 ) nextRandom -= 1.0;
    nextHasBeenSet = true;
  }

  return a;
}


void NonRandomEngine::flatArray(const int size, double* vect) {
  for (int i = 0; i < size; ++i) {
    vect[i] = flat();
  }
}

std::ostream & NonRandomEngine::put (std::ostream & os) const {
  std::string beginMarker = "NonRandomEngine-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
}

std::vector<unsigned long> NonRandomEngine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<NonRandomEngine>());
  std::vector<unsigned long> t;
  v.push_back(static_cast<unsigned long>(nextHasBeenSet));
  v.push_back(static_cast<unsigned long>(sequenceHasBeenSet));
  v.push_back(static_cast<unsigned long>(intervalHasBeenSet));
  t = DoubConv::dto2longs(nextRandom);
  v.push_back(t[0]); v.push_back(t[1]);
  v.push_back(static_cast<unsigned long>(nInSeq));
  t = DoubConv::dto2longs(randomInterval);
  v.push_back(t[0]); v.push_back(t[1]);
  v.push_back(static_cast<unsigned long>(sequence.size()));
  for (unsigned int i=0; i<sequence.size(); ++i) {
    t = DoubConv::dto2longs(sequence[i]);
    v.push_back(t[0]); v.push_back(t[1]);
  }
  return v;
}

std::istream & NonRandomEngine::get (std::istream & is) {
  std::string beginMarker = "NonRandomEngine-begin";
  is >> beginMarker;
  if (beginMarker != "NonRandomEngine-begin") {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\nInput mispositioned or"
	      << "\nNonRandomEngine state description missing or"
	      << "\nwrong engine type found.\n";
    return is;
  }
  return getState(is);
}

std::string NonRandomEngine::beginTag ( )  { 
  return "NonRandomEngine-begin"; 
}  

std::istream & NonRandomEngine::getState (std::istream & is) {
  if ( possibleKeywordInput ( is, "Uvec", nextHasBeenSet ) ) {
    std::vector<unsigned long> v;
    unsigned long uu = 99999;
    unsigned long ssiz = 0;  
    for (unsigned int istart=0; istart < 10; ++istart) {
      is >> uu;
      if (!is) {
	is.clear(std::ios::badbit | is.rdstate());
        std::cout << "istart = " << istart << "\n";
	std::cerr 
	<< "\nNonRandomEngine state (vector) description has no sequence size."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      } 
      v.push_back(uu);
      if (istart==9) ssiz = uu;
    }   
    for (unsigned int ivec=0; ivec < 2*ssiz; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nNonRandomEngine state (vector) description improper."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      }
      v.push_back(uu);
    }
    getState(v);
    return (is);
  }

//  is >> nextHasBeenSet;  Removed, encompassed by possibleKeywordInput()

  std::string  endMarker  = "NonRandomEngine-end";
  is >> sequenceHasBeenSet >> intervalHasBeenSet;
  is >> nextRandom >> nInSeq >> randomInterval;
  unsigned int seqSize;
  is >> seqSize;
  sequence.clear();
  double x;
  for (unsigned int i = 0; i < seqSize; ++i) {
    is >> x;
    sequence.push_back(x);
  }
  is >> endMarker;
  if (endMarker != "NonRandomEngine-end") {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "\n NonRandomEngine state description incomplete."
	      << "\nInput stream is probably mispositioned now." << std::endl;
    return is;
  }
  return is;
}

bool NonRandomEngine::get (const std::vector<unsigned long> & v) {
  if ((v[0] & 0xffffffffUL) != engineIDulong<NonRandomEngine>()) {
    std::cerr << 
    	"\nNonRandomEngine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool NonRandomEngine::getState (const std::vector<unsigned long> & v) {
  unsigned int seqSize = v[9];
  if (v.size() != 2*seqSize + 10 ) {
    std::cerr << 
   "\nNonRandomEngine get:state vector has wrong length - state unchanged\n";
    std::cerr << "  (length = " << v.size() 
              << "; expected " << 2*seqSize + 10 << ")\n"; 
    return false;
  }
  std::vector<unsigned long> t(2);
  nextHasBeenSet     = (v[1]!=0);
  sequenceHasBeenSet = (v[2]!=0);
  intervalHasBeenSet = (v[3]!=0);
  t[0] = v[4]; t[1] = v[5]; nextRandom = DoubConv::longs2double(t);
  nInSeq = v[6];
  t[0] = v[7]; t[1] = v[8]; randomInterval = DoubConv::longs2double(t);
  sequence.clear();
  for (unsigned int i=0; i<seqSize; ++i) {
    t[0] = v[2*i+10]; t[1] = v[2*i+11];
    sequence.push_back(DoubConv::longs2double(t));
  }
  return true;
}


}  // namespace CLHEP

