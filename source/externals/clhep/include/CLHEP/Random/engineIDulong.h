// $Id:$ 
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- engineIDulong ---
//                          function header file
// -----------------------------------------------------------------------

// Class generating new engines from streamed saves.

// =======================================================================
// M Fischler     - Created:  Mar. 8, 2005
// =======================================================================

#ifndef engineIDulong_h
#define engineIDulong_h 1

namespace CLHEP {

unsigned long crc32ul(const std::string & s);

template <class E> 
unsigned long engineIDulong() {
  static const unsigned long id = crc32ul(E::engineName());
  return id;
}

}  // namespace CLHEP

#endif

