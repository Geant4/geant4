// $Id: pyG4HepRep.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4HepRep.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4HepRep.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4HepRep()
{
  class_<G4HepRep, bases<G4VGraphicsSystem> >
    ("G4HepRep", "HepRep visualization module")
    ;
}

