// $Id: pyG4RandomDirection.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RandomDirection.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RandomDirection.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4RandomDirection()
{
  def("G4RandomDirection",  G4RandomDirection);
}

