// $Id: pyG4HepRepFile.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4HepRepFile.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4HepRepFile.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4HepRepFile()
{
  class_<G4HepRepFile, bases<G4VGraphicsSystem> >
    ("G4HepRepFile", "HepRep(File) visualization module")
    ;
}

