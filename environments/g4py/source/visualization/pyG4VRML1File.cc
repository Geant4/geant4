// $Id: pyG4VRML1File.cc,v 1.2 2006-04-25 08:09:46 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VRML1File.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VRML1File.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VRML1File()
{
  class_<G4VRML1File, bases<G4VGraphicsSystem> >
    ("G4VRML1File", "VRML-1(file) visualization module")
    ;
}

