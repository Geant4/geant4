#ifndef CLHEP_MEMORY_H
#define CLHEP_MEMORY_H

// ======================================================================
//
// memory - memory management utilities
//
// ======================================================================

#include <memory>

namespace CLHEP {

template < typename T >
using shared_ptr = std::shared_ptr<T>;
template < typename T >
using weak_ptr = std::weak_ptr<T>;

// ----------------------------------------------------------------------
// do_nothing_deleter - for shared_ptrs not taking ownership
// ----------------------------------------------------------------------

struct do_nothing_deleter {
  inline  void  operator () ( void const * ) const;
};

void
do_nothing_deleter::operator () ( void const * ) const
{ }


}  // namespace CLHEP

#endif  // CLHEP_MEMORY_H
//
// ======================================================================
