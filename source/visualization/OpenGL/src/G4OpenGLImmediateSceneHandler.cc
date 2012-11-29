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
// $Id$
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
#include "G4AttHolder.hh"

#include <typeinfo>

G4OpenGLImmediateSceneHandler::G4OpenGLImmediateSceneHandler
(G4VGraphicsSystem& system,const G4String& name):
  G4OpenGLSceneHandler (system, fSceneIdCount++, name)
{}

G4OpenGLImmediateSceneHandler::~G4OpenGLImmediateSceneHandler ()
{}

#include <iomanip>

G4bool G4OpenGLImmediateSceneHandler::AddPrimitivePreamble(const G4Visible& visible)
{
  const G4Colour& c = GetColour (visible);
  G4double opacity = c.GetAlpha ();
  
  G4bool transparency_enabled = true;
  G4bool isMarkerNotHidden = true;
  G4OpenGLViewer* pViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
  if (pViewer) {
    transparency_enabled = pViewer->transparency_enabled;
    isMarkerNotHidden = pViewer->fVP.IsMarkerNotHidden();
  }
  
  G4bool isMarker = false;
  try {
    (void) dynamic_cast<const G4VMarker&>(visible);
    isMarker = true;
  }
  catch (std::bad_cast) {}
  
  G4bool isPolyline = false;
  try {
    (void) dynamic_cast<const G4Polyline&>(visible);
    isPolyline = true;
  }
  catch (std::bad_cast) {}
  
  G4bool isMarkerOrPolyline = isMarker || isPolyline;
  G4bool treatAsTransparent = transparency_enabled && opacity < 1.;
  G4bool treatAsNotHidden = isMarkerNotHidden && (isMarker || isPolyline);
  
  if (fProcessing2D) glDisable (GL_DEPTH_TEST);
  else {
    if (isMarkerOrPolyline && isMarkerNotHidden)
      glDisable (GL_DEPTH_TEST);
    else {glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LEQUAL);}
  }

  if (fThreePassCapable) {
    
    // Ensure transparent objects are drawn opaque ones and before
    // non-hidden markers.  The problem of blending/transparency/alpha
    // is quite a tricky one - see History of opengl-V07-01-01/2/3.
    if (!(fSecondPassForTransparency || fThirdPassForNonHiddenMarkers)) {
      // First pass...
      if (treatAsTransparent) {  // Request pass for transparent objects...
        fSecondPassForTransparencyRequested = true;
      }
      if (treatAsNotHidden) {    // Request pass for non-hidden markers...
        fThirdPassForNonHiddenMarkersRequested = true;
      }
      // On first pass, transparent objects and non-hidden markers are not drawn...
      if (treatAsTransparent || treatAsNotHidden) {
        return false;
      }
    }
    
    // On second pass, only transparent objects are drawn...
    if (fSecondPassForTransparency) {
      if (!treatAsTransparent) {
        return false;
      }
    }
    
    // On third pass, only non-hidden markers are drawn...
    if (fThirdPassForNonHiddenMarkers) {
      if (!treatAsNotHidden) {
        return false;
      }
    }
  }  // fThreePassCapable
  
  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking()) {
    glLoadName(++fPickName);
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(visible, holder);
    fPickMap[fPickName] = holder;
  }

  if (transparency_enabled) {
    glColor4d(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());
  } else {
    glColor3d(c.GetRed(),c.GetGreen(),c.GetBlue());    
  }

  return true;
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polyline& polyline)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polyline);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(polyline);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polymarker& polymarker)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polymarker);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(polymarker);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Text& text)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Text&).
  G4bool furtherprocessing = AddPrimitivePreamble(text);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(text);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Circle& circle)
{
  G4bool furtherprocessing = AddPrimitivePreamble(circle);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(circle);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Square& square)
{
  G4bool furtherprocessing = AddPrimitivePreamble(square);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(square);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Scale& scale)
{
  G4bool furtherprocessing = AddPrimitivePreamble(scale);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(scale);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4Polyhedron&).
  G4bool furtherprocessing = AddPrimitivePreamble(polyhedron);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(polyhedron);
  }
}

void G4OpenGLImmediateSceneHandler::AddPrimitive (const G4NURBS& nurbs)
{
  // Note: colour is still handled in
  // G4OpenGLSceneHandler::AddPrimitive(const G4NURBS&).
  G4bool furtherprocessing = AddPrimitivePreamble(nurbs);
  if (furtherprocessing) {
    G4OpenGLSceneHandler::AddPrimitive(nurbs);
  }
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

  // See all primitives immediately...  At least soon...
  ScaledFlush();

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
  glDisable(GL_DEPTH_TEST);  // But see parent scene handler!!  In
  glDisable (GL_LIGHTING);   // some cases, we need to re-iterate this.
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

void G4OpenGLImmediateSceneHandler::ClearTransientStore ()
{
  // Nothing to do except redraw the scene ready for the next event.
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

G4int G4OpenGLImmediateSceneHandler::fSceneIdCount = 0;

#endif
