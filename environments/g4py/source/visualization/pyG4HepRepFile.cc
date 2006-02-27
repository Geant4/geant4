// $Id: pyG4HepRepFile.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

