// $Id: pyG4OpenGLStoredX.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4OpenGLStoredX.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_OPENGLX

#include <boost/python.hpp>
#include "G4OpenGLStoredX.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4OpenGLStoredX()
{
  class_<G4OpenGLStoredX, bases<G4VGraphicsSystem> >
    ("G4OpenGLStoredX", "OpenGL(Stored X) visualization module")
    ;
}

#endif
