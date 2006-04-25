// $Id: pyG4OpenGLImmediateX.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4OpenGLImmediateX.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_OPENGLX

#include <boost/python.hpp>
#include "G4OpenGLImmediateX.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4OpenGLImmediateX()
{
  class_<G4OpenGLImmediateX, bases<G4VGraphicsSystem> >
    ("G4OpenGLImmediateX", "OpenGL(Immediate X) visualization module")
    ;
}

#endif

