// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateSceneHandler.cc,v 1.2 1999-01-11 00:47:42 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"

#include "G4OpenGLImmediateSceneHandler.hh"

G4OpenGLImmediateSceneHandler::G4OpenGLImmediateSceneHandler (G4VGraphicsSystem& system,
						const G4String& name):
G4OpenGLSceneHandler (system, fSceneIdCount++, name)
{
  fSceneCount++;
}

G4OpenGLImmediateSceneHandler::~G4OpenGLImmediateSceneHandler ()
{
  fSceneCount--;
}

#include <iomanip.h>

void G4OpenGLImmediateSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  G4VSceneHandler::BeginPrimitives (objectTransformation);
  glPushMatrix();
  G4OpenGLTransform3D oglt (objectTransformation);

  /*************************** Check matrix.
  const GLdouble* m = oglt.GetGLMatrix ();
  G4cout << "G4OpenGLTransform3D matrix:";
  for (int i = 0; i < 16; i++) {
    if ((i % 4) == 0) G4cout << '\n';
    G4cout << setw (15) << m[i];
  }
  G4cout << endl;
  *****************************************/

  glMultMatrixd (oglt.GetGLMatrix ());
}

void G4OpenGLImmediateSceneHandler::EndPrimitives () {
  glPopMatrix();
  if (fReadyForTransients) {
    glFlush ();
  }
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenGLImmediateSceneHandler::BeginModeling () {

  if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
  }
  G4VSceneHandler::BeginModeling();
}

void G4OpenGLImmediateSceneHandler::EndModeling () {
  G4VSceneHandler::EndModeling ();
  if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
    //    glDisable (GL_POLYGON_OFFSET_FILL);
  }
}

G4int G4OpenGLImmediateSceneHandler::GetSceneCount () {
  return fSceneCount;
}

G4int G4OpenGLImmediateSceneHandler::fSceneIdCount = 0;

G4int G4OpenGLImmediateSceneHandler::fSceneCount = 0;

#endif

