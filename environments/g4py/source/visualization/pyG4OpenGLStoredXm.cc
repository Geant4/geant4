// $Id: pyG4OpenGLStoredXm.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4OpenGLStoredXm.cc
//
//                                         2005 Q
// ====================================================================
#ifdef G4VIS_USE_OPENGLXM

#include <boost/python.hpp>
#include "G4OpenGLStoredXm.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4OpenGLStoredXm()
{
  class_<G4OpenGLStoredXm, bases<G4VGraphicsSystem> >
    ("G4OpenGLStoredXm", "OpenGL(Stored XM) visualization module")
    ;
}

#endif
