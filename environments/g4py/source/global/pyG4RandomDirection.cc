// $Id: pyG4RandomDirection.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

