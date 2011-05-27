// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- EngineFactory ---
//                          class header file
// -----------------------------------------------------------------------

// Class generating new engines from streamed saves.

// =======================================================================
// M Fischler     - Created:  12/21/04
// =======================================================================

#ifndef EngineFactory_h
#define EngineFactory_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

class EngineFactory {
public:
  static HepRandomEngine* newEngine(std::istream & is);
  static HepRandomEngine* newEngine(std::vector<unsigned long> const & v);    
};

}  // namespace CLHEP

#endif

