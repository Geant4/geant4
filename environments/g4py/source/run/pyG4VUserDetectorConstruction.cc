// $Id: pyG4VUserDetectorConstruction.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VUserDetectorConstruction.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VUserDetectorConstruction.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VUserDetectorConstruction()
{
  class_<G4VUserDetectorConstruction, boost::noncopyable>
    ("G4VUserDetectorConstruction",
     "base class of user detector construction", no_init)
    ;
}
