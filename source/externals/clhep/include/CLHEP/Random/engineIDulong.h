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
  static G4ThreadLocal unsigned long *id_G4MT_TLS_ = 0 ; if (!id_G4MT_TLS_) {id_G4MT_TLS_ = new  unsigned long  ; *id_G4MT_TLS_= crc32ul(E::engineName()) ; }  unsigned long &id = *id_G4MT_TLS_;
  return id;
}

}  // namespace CLHEP

#endif

