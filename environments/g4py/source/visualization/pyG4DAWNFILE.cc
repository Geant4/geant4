// $Id: pyG4DAWNFILE.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4DAWNFILE.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4DAWNFILE.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4DAWNFILE()
{
  class_<G4DAWNFILE, bases<G4VGraphicsSystem> >
    ("G4DAWNFILE", "DAWN(file) visualization module")
    ;
}

