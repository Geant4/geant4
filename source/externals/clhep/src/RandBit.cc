// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandBit ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// M Fischler     - Created from RandFlat.cc, deleting almost all the content
//		    since inheritance takes care of it.  2/15/00
// M Fischler     - put and get to/from streams 12/10/04
// =======================================================================

#include "CLHEP/Random/RandBit.h"
#include <string>

namespace CLHEP {

std::string RandBit::name() const {return "RandBit";}

RandBit::~RandBit() {
}

std::ostream & RandBit::put ( std::ostream & os ) const {
  os << " " << name() << "\n";
  RandFlat::put(os);
  return os;
}

std::istream & RandBit::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  RandFlat::get(is);
  return is;
}
  
}  // namespace CLHEP

