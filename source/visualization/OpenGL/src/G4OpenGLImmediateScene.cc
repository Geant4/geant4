// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateScene.cc,v 1.1 1999-01-07 16:14:56 gunter Exp $
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

#include "G4OpenGLScene.hh"
#include "G4OpenGLView.hh"
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

#include "G4OpenGLImmediateScene.hh"

G4OpenGLImmediateScene::G4OpenGLImmediateScene (G4VGraphicsSystem& system,
						const G4String& name):
G4OpenGLScene (system, fSceneIdCount++, name)
{
  fSceneCount++;
}

G4OpenGLImmediateScene::~G4OpenGLImmediateScene ()
{
  fSceneCount--;
}

#include <iomanip.h>

void G4OpenGLImmediateScene::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  G4VScene::BeginPrimitives (objectTransformation);
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

void G4OpenGLImmediateScene::EndPrimitives () {
  glPopMatrix();
  if (fReadyForTransients) {
    glFlush ();
  }
  G4VScene::EndPrimitives ();
}

void G4OpenGLImmediateScene::BeginModeling () {

  if (fpView -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
  }
  G4VScene::BeginModeling();
}

void G4OpenGLImmediateScene::EndModeling () {
  G4VScene::EndModeling ();
  if (fpView -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
    //    glDisable (GL_POLYGON_OFFSET_FILL);
  }
}

G4int G4OpenGLImmediateScene::GetSceneCount () {
  return fSceneCount;
}

G4int G4OpenGLImmediateScene::fSceneIdCount = 0;

G4int G4OpenGLImmediateScene::fSceneCount = 0;

#endif

