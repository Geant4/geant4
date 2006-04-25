// $Id: pyG4TouchableHistory.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VTouchable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TouchableHistory.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4TouchableHistory()
{
  class_<G4TouchableHistory, G4TouchableHistory*, bases<G4VTouchable> >
    ("G4TouchableHistory", "touchable history class")
    ;
}

