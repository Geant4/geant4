// $Id: pyG4VRML2File.cc,v 1.2 2006-04-25 08:09:46 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VRML2File.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VRML2File.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VRML2File()
{
  class_<G4VRML2File, bases<G4VGraphicsSystem> >
    ("G4VRML2File", "VRML-2(file) visualization module")
    ;
}

