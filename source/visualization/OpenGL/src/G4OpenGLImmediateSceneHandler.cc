//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLImmediateSceneHandler.cc,v 1.35 2010-06-03 20:35:19 allison Exp $
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

#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4OpenGLViewer.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Scale.hh"
#include "G4Polyhedron.hh"

G4OpenGLImmediateSceneHandler::G4OpenGLImmediateSceneHandler
(G4VGraphicsSystem& system,const G4String& name):
  G4OpenGLSceneHandler (system, fSceneIdCount++, name)
{}

G4OpenGLImmediateSceneHandler::~G4OpenGLImmediateSceneHandler ()
{}

#include <iomanip>

void G4OpenGLImmediateSceneHandler::AddPrimitivePreamble(const G4Visible& visible)
{
  if (fpViewer->GetViewParameters().IsPicking()) {
    glLoadName(++fPickName);
    fPickMap[fPickName] = 0;
  }

  const G4Colour& c = GetColour (visible);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polyline& polyline)
{
  AddPrimitivePreamble(polyline);
  G4OpenGLSceneHandler::AddPrimitive(polyline);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  AddPrimitivePreamble(polymarker);
  G4OpenGLSceneHandler::AddPrimitive(polymarker);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Text& text)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Text&).
  AddPrimitivePreamble(text);
  G4OpenGLSceneHandler::AddPrimitive(text);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Circle& circle)
{
  AddPrimitivePreamble(circle);
  G4OpenGLSceneHandler::AddPrimitive(circle);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Square& square)
{
  AddPrimitivePreamble(square);
  G4OpenGLSceneHandler::AddPrimitive(square);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Scale& scale)
{
  AddPrimitivePreamble(scale);
  G4OpenGLSceneHandler::AddPrimitive(scale);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Polyhedron&).
  AddPrimitivePreamble(polyhedron);
  G4OpenGLSceneHandler::AddPrimitive(polyhedron);
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4NURBS& nurbs)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4NURBS&).
  AddPrimitivePreamble(nurbs);
  G4OpenGLSceneHandler::AddPrimitive(nurbs);
}

void G4OpenGLImmediateSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation)
{
  G4OpenGLSceneHandler::BeginPrimitives (objectTransformation);

  G4OpenGLTransform3D oglt (objectTransformation);

  glPushMatrix();

  /*************************** Check matrix.
  const GLdouble* m = oglt.GetGLMatrix ();
  G4cout << "G4OpenGLTransform3D matrix:";
  for (int i = 0; i < 16; i++) {
    if ((i % 4) == 0) G4cout << '\n';
    G4cout << std::setw (15) << m[i];
  }
  G4cout << G4endl;
  *****************************************/

  glMultMatrixd (oglt.GetGLMatrix ());
}

void G4OpenGLImmediateSceneHandler::EndPrimitives ()
{
  glPopMatrix();

  // See all primitives immediately...
  glFlush ();

  G4OpenGLSceneHandler::EndPrimitives ();
}

void G4OpenGLImmediateSceneHandler::BeginPrimitives2D
(const G4Transform3D& objectTransformation)
{
  G4OpenGLSceneHandler::BeginPrimitives2D(objectTransformation);

  // Push current 3D world matrices and load identity to define screen
  // coordinates...
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho (-1., 1., -1., 1., -G4OPENGL_FLT_BIG, G4OPENGL_FLT_BIG);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  G4OpenGLTransform3D oglt (objectTransformation);
  glMultMatrixd (oglt.GetGLMatrix ());
}

void G4OpenGLImmediateSceneHandler::EndPrimitives2D()
{
  // Pop current 3D world matrices back again...
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();

  // See all primitives immediately...
  glFlush ();

  G4OpenGLSceneHandler::EndPrimitives2D ();
}

void G4OpenGLImmediateSceneHandler::BeginModeling () {
  G4VSceneHandler::BeginModeling();
}

void G4OpenGLImmediateSceneHandler::EndModeling () {
  G4VSceneHandler::EndModeling ();
}

void G4OpenGLImmediateSceneHandler::ClearTransientStore () {

  G4VSceneHandler::ClearTransientStore ();

  // Make sure screen corresponds to graphical database...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

void G4OpenGLImmediateSceneHandler::RequestPrimitives (const G4VSolid& solid)
{
  if (fReadyForTransients) {
    // Always draw transient solids, e.g., hits represented as solids.
    // (As we have no control over the order of drawing of transient
    // objects, we cannot do anything about transparent ones, as
    // below, so always draw them.)
    G4VSceneHandler::RequestPrimitives (solid);
    return;
  }

  // For non-transient (run-duration) objects, ensure transparent
  // objects are drawn last.  The problem of
  // blending/transparency/alpha is quite a tricky one - see History
  // of opengl-V07-01-01/2/3.
  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes(fpVisAttribs);
  const G4Colour& c = pVA -> GetColour ();
  G4double opacity = c.GetAlpha ();

  if (!fSecondPass) {
    G4bool transparency_enabled = true;
    G4OpenGLViewer* pViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
    if (pViewer) transparency_enabled = pViewer->transparency_enabled;
    if (transparency_enabled && opacity < 1.) {
      // On first pass, transparent objects are not drawn, but flag is set...
      fSecondPassRequested = true;
      return;
    }
  }

  // On second pass, opaque objects are not drwan...
  if (fSecondPass && opacity >= 1.) return;

  // Else invoke base class method...
  G4VSceneHandler::RequestPrimitives (solid);
}

G4int G4OpenGLImmediateSceneHandler::fSceneIdCount = 0;

#endif
