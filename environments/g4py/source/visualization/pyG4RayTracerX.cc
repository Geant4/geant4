// $Id: pyG4RayTracerX.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RayTracerX.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_RAYTRACERX

#include "pyG4Version.hh"

#if G4VERSION_NUMBER >= 800

#include <boost/python.hpp>
#include "G4RayTracerX.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4RayTracerX()
{
  class_<G4RayTracerX, bases<G4VGraphicsSystem> >
    ("G4RayTracerX", "RayTracerX visualization module")
    ;
}

#endif
#endif
