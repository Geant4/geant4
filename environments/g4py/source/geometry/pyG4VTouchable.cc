// $Id: pyG4VTouchable.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VTouchable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VTouchable.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VTouchable()
{
  class_<G4VTouchable, G4VTouchable*, boost::noncopyable>
    ("G4VTouchable", "touchable class", no_init)
    ;
}

