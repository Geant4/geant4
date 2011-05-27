// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                      --- StaticRandomStates ---
//                      class implementation file
// -----------------------------------------------------------------------
//
// =======================================================================
// Mark Fischler  - Created: Dec. 21, 2004
// Mark Fischler  - Modified restore() to utilize anonymous engine input
//                  to create anonymous restore of the static distributions
//
// =======================================================================

#include "CLHEP/Random/StaticRandomStates.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include <string>
#include <sstream>

//======================//
//                      //
// Maintenance warning: //
//			//
//======================//
//
// Currently, only two distributions (RandFlat and RandGauss) have cached
// distribution state.  All such distributions must be saved below, so if
// another such distribution is added, this implementation file must be 
// modified to reflect that.

namespace CLHEP {


std::ostream & StaticRandomStates::save(std::ostream & os){
  RandGauss::saveFullState(os);
  RandFlat::saveDistState(os);
  return os;
}

#ifdef NOTYET
std::istream & StaticRandomStates::restore(std::istream & is) {
  RandGauss::restoreFullState(is);
  RandFlat::restoreDistState(is);
  return is;
}
#endif

std::istream & StaticRandomStates::restore(std::istream & is) {
  HepRandomEngine * e = HepRandom::getTheEngine();
  HepRandomEngine *ne = HepRandomEngine::newEngine(is);
  if ( !is ) return is;
  if ( !ne ) return is;
  if (ne->name() == e->name()) {
    // Because e has const data members, cannot simply do *e = *ne
    std::ostringstream os;
    os << *ne;
    std::istringstream istst(os.str());
    istst >> *e;
    if (!istst) {
      std::cerr << "???? Unexpected behavior in StaticRandomStates::restore:\n"
        << "The new engine, which had been input successfully from istream\n"
	<< "has encountered a problem when used to set state of theEngine\n";
      is.clear(std::ios::badbit | is.rdstate());
      return is;
    }
  } else {
    HepRandom::setTheEngine(ne);
  }
  RandGauss::restoreDistState(is);
  RandFlat::restoreDistState(is);
  return is;
}

}  // namespace CLHEP
