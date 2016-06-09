//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLSceneHandler.cc,v 1.32 2005/04/17 16:08:43 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL stored scene - creates OpenGL display lists.
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLFontBaseStore.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4VMarker.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

G4OpenGLSceneHandler::G4OpenGLSceneHandler (G4VGraphicsSystem& system,
			      G4int id,
			      const G4String& name):
G4VSceneHandler (system, id, name)
{}

G4OpenGLSceneHandler::~G4OpenGLSceneHandler ()
{
  ClearStore ();
}

const GLubyte G4OpenGLSceneHandler::fStippleMaskHashed [128] = {
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55,
  0x55,0x55,0x55,0x55,0x55,0x55,0x55,0x55
};

//Method for handling G4Polyline objects (from tracking or wireframe).
void G4OpenGLSceneHandler::AddPrimitive (const G4Polyline& line)
{
  G4int nPoints = line.size ();
  if (nPoints <= 0) return;

  const G4Colour& c = GetColour (line);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());

  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else glEnable (GL_DEPTH_TEST);

  glDisable (GL_LIGHTING);
  glBegin (GL_LINE_STRIP);

  for (G4int iPoint = 0; iPoint < nPoints; iPoint++) {
  G4double x, y, z;
    x = line[iPoint].x(); 
    y = line[iPoint].y();
    z = line[iPoint].z();
    glVertex3d (x, y, z);
  }
  glEnd ();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) {

  const G4Colour& c = GetColour (text);  // Picks up default if none.
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (text, sizeType);
  G4ThreeVector position (*fpObjectTransformation * text.GetPosition ());
  G4String textString = text.GetText();

  G4int font_base = G4OpenGLFontBaseStore::GetFontBase(fpViewer,size);
  if (font_base < 0) {
    static G4int callCount = 0;
    ++callCount;
    if (callCount <= 10 || callCount%100 == 0) {
      G4cout <<
	"G4OpenGLSceneHandler::AddPrimitive (const G4Text&) call count "
	     << callCount <<
	"\n  No fonts available."
	"\n  Called with text \""
	     << text.GetText ()
	     << "\"\n  at " << position
	     << ", size " << size
	     << ", offsets " << text.GetXOffset () << ", " << text.GetYOffset ()
	     << ", type " << G4int(sizeType)
	     << ", colour " << c
	     << G4endl;
    }
    return;
  }
  const char* textCString = textString.c_str();
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  glDisable (GL_DEPTH_TEST);
  glDisable (GL_LIGHTING);
  
  glRasterPos3f(position.x(),position.y(),position.z());
  // No action on offset or layout at present.
  glPushAttrib(GL_LIST_BIT);
  glListBase(font_base);
  glCallLists(strlen(textCString), GL_UNSIGNED_BYTE, (GLubyte *)textCString);
  glPopAttrib();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Circle& circle) {
  glEnable (GL_POINT_SMOOTH);
  AddCircleSquare (circle, 24);
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Square& square) {
  glDisable (GL_POINT_SMOOTH);
  AddCircleSquare (square, 4);
}

void G4OpenGLSceneHandler::AddCircleSquare
(const G4VMarker& marker,
 G4int nSides) {

  const G4Colour& c = GetColour (marker);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  
  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else glEnable (GL_DEPTH_TEST);
  
  glDisable (GL_LIGHTING);
  
  G4VMarker::FillStyle style = marker.GetFillStyle();
  
  switch (style) {
  case G4VMarker::noFill: 
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    break;
    
  case G4VMarker::hashed:
    /*
    G4cout << "Hashed fill style in G4OpenGLSceneHandler."
	   << "\n  Not implemented.  Using G4VMarker::filled."
	   << G4endl;
    */
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glPolygonStipple (fStippleMaskHashed);
    // See also:
    //   if (style == G4VMarker::filled || style == G4VMarker::hashed)...
    // (twice) below.
    break;
    
  case G4VMarker::filled:
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;
    
  default:
    G4cout << "Unrecognised fill style in G4OpenGLSceneHandler."
	   << "\n  Using G4VMarker::filled."
	   << G4endl;
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;
    
  }

  // A few useful quantities...
  G4Point3D centre = marker.GetPosition();
  G4bool userSpecified = (marker.GetWorldSize() || marker.GetScreenSize());
  const G4VMarker& def = fpViewer -> GetViewParameters().GetDefaultMarker();
  const G4Vector3D& viewpointDirection =
    fpViewer -> GetViewParameters().GetViewpointDirection();
  const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  G4double scale = fpViewer -> GetViewParameters().GetGlobalMarkerScale();
  G4double size = scale *
    userSpecified ? marker.GetWorldSize() : def.GetWorldSize();

  // Find "size" of marker in world space (but see note below)...
  G4double worldSize;
  if (size) {  // Size specified in world coordinates.
    worldSize = size;
  }
  else { // Size specified in screen (window) coordinates.

    // Find window coordinates of centre...
    GLdouble modelMatrix[16];
    glGetDoublev (GL_MODELVIEW_MATRIX, modelMatrix);
    G4double projectionMatrix[16];
    glGetDoublev (GL_PROJECTION_MATRIX, projectionMatrix);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    GLdouble winx, winy, winz;
    gluProject(centre.x(), centre.y(), centre.z(),
	       modelMatrix, projectionMatrix, viewport,
	       &winx, &winy, &winz);

    // Determine ratio window:world...
    const G4Vector3D inScreen = (up.cross(viewpointDirection)).unit();
    const G4Vector3D p = centre + inScreen;
    GLdouble winDx, winDy, winDz;
    gluProject(p.x(), p.y(), p.z(),
               modelMatrix, projectionMatrix, viewport,
               &winDx, &winDy, &winDz);
    G4double winWorldRatio = std::sqrt(std::pow(winx - winDx, 2) +
				  std::pow(winy - winDy, 2));
    G4double winSize = scale *
      userSpecified ? marker.GetScreenSize() : def.GetScreenSize();
    worldSize = winSize / winWorldRatio;
  }

  // Draw...
  DrawXYPolygon (worldSize, centre, nSides);
}

/***************************************************
Note: We have to do it this way round so that when a global
transformation is applied, such as with /vis/viewer/set/viewpoint,
the markers follow the world coordinates without having to
recreate the display lists.  The down side is that the markers
rotate.  The only way to avoid this is to play with the modelview
and projection matrices of OpenGL - which I need to think about.
For future reference, here is the code to draw in window
coordinates; its down side is that markers do not follow global
transformations.  Some clever stuff is needed.

  ...
  // Find window coordinates of centre...
  GLdouble modelMatrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, modelMatrix);
  G4double projectionMatrix[16];
  glGetDoublev (GL_PROJECTION_MATRIX, projectionMatrix);
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);
  GLdouble winx, winy, winz;
  gluProject(centre.x(), centre.y(), centre.z(),
             modelMatrix, projectionMatrix, viewport,
             &winx, &winy, &winz);

  // Find window size...
  G4double winSize;
  if (size) {  // Size specified in world coordinates.
    // Determine size in window coordinates...
    (Note: improve this by using an inScreen vector as above.)
    GLdouble winx1, winy1, winz1;
    gluProject(centre.x() + size, centre.y() + size, centre.z() + size,
               modelMatrix, projectionMatrix, viewport,
               &winx1, &winy1, &winz1);
    winSize = std::sqrt((std::pow(winx - winx1, 2) +
                    std::pow(winy - winy1, 2) +
                    std::pow(winz - winz1, 2)) / 3.);
  }
  else {
    winSize = scale *
      userSpecified ? marker.GetScreenSize() : def.GetScreenSize();
  }

  // Prepare to draw in window coordinates...
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(GLdouble(viewport[0]),
             GLdouble(viewport[0] + viewport[2]),
             GLdouble(viewport[1]),
             GLdouble(viewport[1] + viewport[3]));
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Draw in window coordinates...
  DrawScreenPolygon (winSize, G4Point3D(winx, winy, winz), nSides);

  // Re-instate matrices...
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  ...
}

void G4OpenGLSceneHandler::DrawScreenPolygon
(G4double size,
 const G4Point3D& centre,
 G4int nSides) {
  glBegin (GL_POLYGON);
  const G4double dPhi = twopi / nSides;
  const G4double r = size / 2.;
  G4double phi;
  G4int i;
  for (i = 0, phi = -dPhi / 2.; i < nSides; i++, phi += dPhi) {
    G4double x, y, z;
    x = centre.x() + r * std::cos(phi);
    y = centre.y() + r * std::sin(phi);
    z = centre.z();
    glVertex3d (x, y, z);
  }
  glEnd ();
}
**********************************************/

void G4OpenGLSceneHandler::DrawXYPolygon
(G4double size,
 const G4Point3D& centre,
 G4int nSides) {
  const G4Vector3D& viewpointDirection =
    fpViewer -> GetViewParameters().GetViewpointDirection();
  const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  const G4double dPhi = twopi / nSides;
  const G4double radius = size / 2.;
  G4Vector3D start = radius * (up.cross(viewpointDirection)).unit();
  G4double phi;
  G4int i;
  glBegin (GL_POLYGON);
  for (i = 0, phi = -dPhi / 2.; i < nSides; i++, phi += dPhi) {
    G4Vector3D r = start; r.rotate(phi, viewpointDirection);
    G4Vector3D p = centre + r;
    glVertex3d (p.x(), p.y(), p.z());
  }
  glEnd ();
}

//Method for handling G4Polyhedron objects for drawing solids.
void G4OpenGLSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

  //Assume all facets are convex quadrilaterals.
  //Draw each G4Facet individually
  
  if (polyhedron.GetNoFacets() == 0) return;

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (polyhedron.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);
  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);
  
  //Get colour, etc..
  const G4Colour& c = pVA -> GetColour ();
  GLfloat materialColour [4];
  materialColour [0] = c.GetRed ();
  materialColour [1] = c.GetGreen ();
  materialColour [2] = c.GetBlue ();
  materialColour [3] = 1.0;
  GLdouble clear_colour[4];
  glGetDoublev (GL_COLOR_CLEAR_VALUE, clear_colour);

  switch (drawing_style) {
  case (G4ViewParameters::hlhsr):
    // Set up as for hidden line removal but paint polygons...
  case (G4ViewParameters::hlr):
    glEnable (GL_STENCIL_TEST);
    // The stencil buffer is cleared in G4OpenGLViewer::ClearView.
    // The procedure below leaves it clear.
    glStencilFunc (GL_ALWAYS, 0, 1);
    glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);    
    glEnable (GL_CULL_FACE);
    glCullFace (GL_BACK);
    glPolygonMode (GL_FRONT, GL_LINE);
    glDisable (GL_LIGHTING);
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    break;
  case (G4ViewParameters::hsr):
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_LESS);    
    glEnable (GL_CULL_FACE);
    glCullFace (GL_BACK); 
    glEnable (GL_LIGHTING);
    glPolygonMode (GL_FRONT, GL_FILL);
    glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
    break;
  case (G4ViewParameters::wireframe):
  default:
    //glDisable (GL_DEPTH_TEST);
    glEnable (GL_DEPTH_TEST);
    glDepthFunc (GL_ALWAYS);    
    glDisable (GL_CULL_FACE);
    glDisable (GL_LIGHTING);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    break;
  }

  //Loop through all the facets...
  glBegin (GL_QUADS);
  G4bool notLastFace;
  do {

    //First, find surface normal and note "not last facet"...
    G4Normal3D SurfaceUnitNormal;
    notLastFace = polyhedron.GetNextUnitNormal (SurfaceUnitNormal);

    //Loop through the four edges of each G4Facet...
    G4bool notLastEdge;
    G4Point3D vertex[4];
    G4int edgeFlag[4];
    G4int edgeCount = 0;
    glNormal3d(SurfaceUnitNormal.x(), SurfaceUnitNormal.y(),
	       SurfaceUnitNormal.z());
    do {
      notLastEdge = polyhedron.GetNextVertex (vertex[edgeCount], 
					      edgeFlag[edgeCount]);
      // Check to see if edge is visible or not...
      if (isAuxEdgeVisible) {
	edgeFlag[edgeCount] = G4int (true);
      }
      if (edgeFlag[edgeCount]) {
	glEdgeFlag (GL_TRUE);
      } else {
	glEdgeFlag (GL_FALSE);
      }
      glVertex3d (vertex[edgeCount].x(), 
		  vertex[edgeCount].y(),
		  vertex[edgeCount].z());
      edgeCount++;
    } while (notLastEdge && edgeCount < 4);
    // HEPPolyhedron produces triangles too; in that case add an extra vertex..
    while (edgeCount < 4) {
      vertex[edgeCount] = vertex[edgeCount-1];
      edgeFlag[edgeCount] = G4int (false);
      glEdgeFlag (GL_FALSE);
      glVertex3d (vertex[edgeCount].x(),
		  vertex[edgeCount].y(), 
		  vertex[edgeCount].z());
      edgeCount++;
    }
    // Trap situation where number of edges is > 4...
    while (notLastEdge) {
      notLastFace = polyhedron.GetNextUnitNormal (SurfaceUnitNormal);
      edgeCount++;
    }
    if (edgeCount > 4) {
      G4cerr <<
	"G4OpenGLSceneHandler::AddPrimitive(G4Polyhedron): WARNING";
      if (fpCurrentPV) {
	G4cerr <<
	"\n  Volume " << fpCurrentPV->GetName() <<
	", Solid " << fpCurrentLV->GetSolid()->GetName() <<
	  " (" << fpCurrentLV->GetSolid()->GetEntityType();
      }
      G4cerr<<
	"\n   G4Polyhedron facet with " << edgeCount << " edges" << G4endl;
    }

    // Do it all over again (twice) for hlr...
    if  (drawing_style == G4ViewParameters::hlr ||
	 drawing_style == G4ViewParameters::hlhsr) {

      glEnd ();  // Placed here to balance glBegin above, allowing GL
		 // state changes below, then glBegin again.  Avoids
		 // having glBegin/End pairs *inside* loop in the more
		 // usual case of no hidden line removal.

      // Draw through stencil...
      glStencilFunc (GL_EQUAL, 0, 1);
      glStencilOp (GL_KEEP, GL_KEEP, GL_KEEP);
      glPolygonMode (GL_FRONT, GL_FILL);
      if (drawing_style == G4ViewParameters::hlhsr) glEnable (GL_LIGHTING);
      glBegin (GL_QUADS);
      if  (drawing_style == G4ViewParameters::hlr) {
	glColor4dv (clear_colour);
      } else {  // drawing_style == G4ViewParameters::hlhsr
	glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
      }
      glNormal3d(SurfaceUnitNormal.x(), SurfaceUnitNormal.y(),
		 SurfaceUnitNormal.z());
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
	if (edgeFlag[edgeCount]) {
	  glEdgeFlag (GL_TRUE);
	} else {
	  glEdgeFlag (GL_FALSE);
	}
	glVertex3d (vertex[edgeCount].x(), 
		    vertex[edgeCount].y(),
		    vertex[edgeCount].z());
      }
      glEnd ();
      glDisable (GL_LIGHTING);

      // and once more to reset the stencil bits...
      glStencilFunc (GL_ALWAYS, 0, 1);
      glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
      glDepthFunc (GL_LEQUAL);  // to make sure line gets drawn.  
      glPolygonMode (GL_FRONT, GL_LINE);
      glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
      glBegin (GL_QUADS);
      for (int edgeCount = 0; edgeCount < 4; ++edgeCount) {
	if (edgeFlag[edgeCount]) {
	  glEdgeFlag (GL_TRUE);
	} else {
	  glEdgeFlag (GL_FALSE);
	}
	glVertex3d (vertex[edgeCount].x(), 
		    vertex[edgeCount].y(),
		    vertex[edgeCount].z());
      }
      glEnd ();
      glDepthFunc (GL_LESS);   // Revert for next quadrilateral.
      glBegin (GL_QUADS);      // Ready for next quadrilateral.  GL
			       // says it ignores incomplete
			       // quadrilaterals, so final empty
			       // glBegin/End sequence should be OK.
    }
  } while (notLastFace);  
  
  glEnd ();
  glDisable (GL_STENCIL_TEST);
}

//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4OpenGLSceneHandler::AddPrimitive (const G4NURBS& nurb) {

  GLUnurbsObj *gl_nurb;
  gl_nurb = gluNewNurbsRenderer ();

  GLfloat *u_knot_array, *u_knot_array_ptr;
  u_knot_array = u_knot_array_ptr = new GLfloat [nurb.GetnbrKnots(G4NURBS::U)];
  G4NURBS::KnotsIterator u_iterator (nurb, G4NURBS::U);
  while (u_iterator.pick (u_knot_array_ptr++));

  GLfloat *v_knot_array, *v_knot_array_ptr;
  v_knot_array = v_knot_array_ptr = new GLfloat [nurb.GetnbrKnots(G4NURBS::V)];
  G4NURBS::KnotsIterator v_iterator (nurb, G4NURBS::V);
  while (v_iterator.pick (v_knot_array_ptr++));

  GLfloat *ctrl_pnt_array, *ctrl_pnt_array_ptr;
  ctrl_pnt_array = ctrl_pnt_array_ptr =
    new GLfloat [nurb.GettotalnbrCtrlPts () * G4NURBS::NofC];
  G4NURBS::CtrlPtsCoordsIterator c_p_iterator (nurb);
  while (c_p_iterator.pick (ctrl_pnt_array_ptr++));

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (nurb.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);
  //G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);
  
  //Get colour, etc..
  const G4Colour& c = pVA -> GetColour ();

  switch (drawing_style) {

  case (G4ViewParameters::hlhsr):
    //    G4cout << "Hidden line removal not implememented in G4OpenGL.\n"
    // << "Using hidden surface removal." << G4endl;
  case (G4ViewParameters::hsr):
    {
      glEnable (GL_LIGHTING);
      glEnable (GL_DEPTH_TEST);
      glEnable (GL_AUTO_NORMAL);
      glEnable (GL_NORMALIZE);
      gluNurbsProperty (gl_nurb, GLU_DISPLAY_MODE, GLU_FILL);
      gluNurbsProperty (gl_nurb, GLU_SAMPLING_TOLERANCE, 50.0);
      GLfloat materialColour [4];
      materialColour [0] = c.GetRed ();
      materialColour [1] = c.GetGreen ();
      materialColour [2] = c.GetBlue ();
      materialColour [3] = 1.0;
      glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
      break;
    }
  case (G4ViewParameters::hlr):
    //    G4cout << "Hidden line removal not implememented in G4OpenGL.\n"
    // << "Using wireframe." << G4endl;
  case (G4ViewParameters::wireframe):
  default:
    glDisable (GL_LIGHTING);
//    glDisable (GL_DEPTH_TEST);
    glEnable (GL_DEPTH_TEST);
    glDisable (GL_AUTO_NORMAL);
    glDisable (GL_NORMALIZE);
    gluNurbsProperty (gl_nurb, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);
    gluNurbsProperty (gl_nurb, GLU_SAMPLING_TOLERANCE, 50.0);
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    break;
  }	

  gluBeginSurface (gl_nurb);
  G4int u_stride = 4;
  G4int v_stride = nurb.GetnbrCtrlPts(G4NURBS::U) * 4;

  gluNurbsSurface (gl_nurb, 
		   nurb.GetnbrKnots (G4NURBS::U), (GLfloat*)u_knot_array, 
		   nurb.GetnbrKnots (G4NURBS::V), (GLfloat*)v_knot_array,
		   u_stride,
		   v_stride,  
		   ctrl_pnt_array,
		   nurb.GetUorder (),
		   nurb.GetVorder (), 
		   GL_MAP2_VERTEX_4);
  
  gluEndSurface (gl_nurb);

  delete [] u_knot_array;  // These should be allocated with smart allocators
  delete [] v_knot_array;  // to avoid memory explosion.
  delete [] ctrl_pnt_array;

  gluDeleteNurbsRenderer (gl_nurb);
}

void G4OpenGLSceneHandler::AddCompound(const G4VTrajectory& traj) {
  G4VSceneHandler::AddCompound(traj);  // For now.
}

void G4OpenGLSceneHandler::AddCompound(const G4VHit& hit) {
  G4VSceneHandler::AddCompound(hit);  // For now.
}

#endif
