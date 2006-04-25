// $Id: pyG4RayTracer.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RayTracer.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RayTracer.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4RayTracer()
{
  class_<G4RayTracer, bases<G4VGraphicsSystem> >
    ("G4RayTracer", "RayTracer visualization module")
    ;
}

