// $Id: pyG4ASCIITree.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ASCIITree.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ASCIITree.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ASCIITree()
{
  class_<G4ASCIITree, bases<G4VGraphicsSystem> >
    ("G4ASCIITree", "ASCII tree visualization module")
    ;
}

