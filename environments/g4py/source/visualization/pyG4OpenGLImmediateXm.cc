// $Id: pyG4OpenGLImmediateXm.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4OpenGLImmediateXm.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_OPENGLXM

#include <boost/python.hpp>
#include "G4OpenGLImmediateXm.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4OpenGLImmediateXm()
{
  class_<G4OpenGLImmediateXm, bases<G4VGraphicsSystem> >
    ("G4OpenGLImmediateXm", "OpenGL(Immediate XM) visualization module")
    ;
}

#endif

