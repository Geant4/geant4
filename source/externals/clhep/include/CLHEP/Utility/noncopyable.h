#ifndef CLHEP_NONCOPYABLE_H
#define CLHEP_NONCOPYABLE_H

// ======================================================================
//
// noncopyable - classes directly/indirectly inheriting won't be copyable
//
// Author:  W. E. Brown; 2010-03-05
//
// ======================================================================


#include "CLHEP/Utility/defs.h"


namespace CLHEP {

class noncopyable
{
protected:
  noncopyable () throw () { }
  ~noncopyable() throw () { }

private:
  noncopyable              ( noncopyable const & );  // = delete;
  noncopyable & operator = ( noncopyable const & );  // = delete;
};  // noncopyable

}  // namespace CLHEP

#endif  // HEP_NONCOPYABLE_H
//
// ======================================================================
