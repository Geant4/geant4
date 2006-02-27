// $Id: pyG4Field.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Field.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Field.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Field()
{
  class_<G4Field, G4Field*, boost::noncopyable>
    ("G4Field", "base class of field", no_init)
    ;
}

