// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLSceneHandler.cc,v 1.6 2001-01-16 18:29:58 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4VMarker.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"

G4OpenGLSceneHandler::G4OpenGLSceneHandler (G4VGraphicsSystem& system,
			      G4int id,
			      const G4String& name):
G4VSceneHandler (system, id, name)
{
  initialize_hlr = true;
  //  glPolygonOffset (1.0, 2);
}

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
  const G4Colour& c = GetColour (line);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());

  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else glEnable (GL_DEPTH_TEST);

  glDisable (GL_LIGHTING);
  glBegin (GL_LINE_STRIP);

  G4int nPoints = line.entries ();
  for (G4int iPoint = 0; iPoint < nPoints; iPoint++) {
  G4double x, y, z;
    x = line(iPoint).x(); 
    y = line(iPoint).y();
    z = line(iPoint).z();
    glVertex3d (x, y, z);
  }
  glEnd ();
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) {
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (text, sizeType);
  G4cout
    << "G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) not implemented yet."
    << "\n  Called with text \"" << text.GetText ()
    << "\" at " << text.GetPosition ()
    << ", size " << size
    << ", offsets " << text.GetXOffset () << ", " << text.GetYOffset ()
    << ", type " << sizeType
    << G4endl;
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
    GLdouble* modelMatrix = new G4double[16];
    glGetDoublev (GL_MODELVIEW_MATRIX, modelMatrix);
    G4double* projectionMatrix = new G4double[16];
    glGetDoublev (GL_PROJECTION_MATRIX, projectionMatrix);
    GLint* viewport = new G4int[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    GLdouble winx, winy, winz;
    gluProject(centre.x(), centre.y(), centre.z(),
	       modelMatrix, projectionMatrix, viewport,
	       &winx, &winy, &winz);

    // Determine ratio window:world...
    const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
    const G4Vector3D inScreen = (up.cross(viewpointDirection)).unit();
    const G4Vector3D p = centre + inScreen;
    GLdouble winDx, winDy, winDz;
    gluProject(p.x(), p.y(), p.z(),
               modelMatrix, projectionMatrix, viewport,
               &winDx, &winDy, &winDz);
    G4double winWorldRatio = sqrt((pow(winx - winDx, 2) +
				   pow(winy - winDy, 2)) / 2.);
    G4double winSize = scale *
      userSpecified ? marker.GetScreenSize() : def.GetScreenSize();
    worldSize = winSize / winWorldRatio;

    delete[] viewport;
    delete[] projectionMatrix;
    delete[] modelMatrix;
  }

  // Draw...
  DrawXYPolygon (worldSize,
		 G4Point3D(centre.x(), centre.y(), centre.z()),
		 viewpointDirection,
		 nSides);
}

/***************************************************
Note: We have to do it this way round so that when a global
transformation is applied, such as with /vis/viewer/set/viewpoint,
the markers follow the world coordinates without having to
recreate the display lists.  The down side is that the markers
rotate.  The only way to avoid this is to play with the modelview
and projection matrices of OpenGL - which I need to think about.
For future reference, here is the code to draw in window
coordinates; it's down side is tha markers do not follow global
transformations.  Some clever stuff is needed.

  ...
  // Find window coordinates of centre...
  GLdouble* modelMatrix = new G4double[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, modelMatrix);
  G4double* projectionMatrix = new G4double[16];
  glGetDoublev (GL_PROJECTION_MATRIX, projectionMatrix);
  GLint* viewport = new G4int[4];
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
    winSize = sqrt((pow(winx - winx1, 2) +
                    pow(winy - winy1, 2) +
                    pow(winz - winz1, 2)) / 3.);
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
  DrawXYPolygon (winSize, G4Point3D(winx, winy, winz), nSides);

  // Re-instate matrices...
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();

  delete[] viewport;
  delete[] projectionMatrix;
  delete[] modelMatrix;
  ...

void G4OpenGLSceneHandler::DrawXYPolygon
(G4double size,
 const G4Point3D& centre,
 G4int nSides) {
  glBegin (GL_POLYGON);
  const G4double dPhi = 2. * M_PI / nSides;
  const G4double r = size / 2.;
  G4double phi;
  G4int i;
  for (i = 0, phi = -dPhi / 2.; i < nSides; i++, phi += dPhi) {
    G4double x, y, z;
    x = centre.x() + r * cos(phi);
    y = centre.y() + r * sin(phi);
    z = centre.z();
    glVertex3d (x, y, z);
  }
  glEnd ();
}
**********************************************/

void G4OpenGLSceneHandler::DrawXYPolygon
(G4double size,
 const G4Point3D& centre,
 const G4Vector3D& normal,
 G4int nSides) {
  const G4Vector3D& up = fpViewer->GetViewParameters().GetUpVector();
  const G4double dPhi = 2. * M_PI / nSides;
  const G4double radius = size / 2.;
  G4Vector3D start = radius * (up.cross(normal)).unit();
  G4double phi;
  G4int i;
  glBegin (GL_POLYGON);
  for (i = 0, phi = -dPhi / 2.; i < nSides; i++, phi += dPhi) {
    G4Vector3D r = start; r.rotate(phi, normal);
    G4Vector3D p = centre + r;
    glVertex3d (p.x(), p.y(), p.z());
  }
  glEnd ();
}

//Method for handling G4Polyhedron objects for drawing solids.
void G4OpenGLSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

  //Assume all facets are convex quadrilaterals.
  //Draw each G4Facet individually
  
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (polyhedron);
  
  //Get colour, etc..
  const G4Colour& c = GetColour (polyhedron);
  
  switch (drawing_style) {
    
    //  case (G4ViewParameters::hlhsr):
    //    cout << "Hidden edge not implememented in G4OpenGL.\n"
    // << "Using hidden surface removal." << G4endl;
    
  case (G4ViewParameters::hsr):
    {
      glPolygonMode (GL_FRONT, GL_FILL);
      glCullFace (GL_BACK); 
      //glEnable (GL_CULL_FACE);
      glEnable (GL_LIGHTING);
      glEnable (GL_DEPTH_TEST);
      GLfloat materialColour [4];
      materialColour [0] = c.GetRed ();
      materialColour [1] = c.GetGreen ();
      materialColour [2] = c.GetBlue ();
      materialColour [3] = 1.0;
      glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColour);
      break;
    }
    
  case (G4ViewParameters::hlr):
    if (initialize_hlr) {
      glEnable (GL_STENCIL_TEST);
      glClear (GL_STENCIL_BUFFER_BIT);
      glStencilFunc (GL_ALWAYS, 0, 1);
      glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
      glGetDoublev (GL_COLOR_CLEAR_VALUE, clear_colour);
      initialize_hlr = false;
      //      glEnable (GL_POLYGON_OFFSET_FILL);
    }
      glDepthFunc (GL_LESS);    
    //drop through to wireframe for first pass...
    
  case (G4ViewParameters::wireframe):
    
  default:
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glDisable (GL_CULL_FACE);
    glDisable (GL_LIGHTING);
    glEnable (GL_DEPTH_TEST);
    glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
    break;
  }	
  
  glBegin (GL_QUADS);
  
  //Loop through all the facets...
  G4bool notLastFace;
  G4Normal3D SurfaceUnitNormal;
  do {
    
    //First, find surface normal for the facet...
    notLastFace = polyhedron.GetNextUnitNormal (SurfaceUnitNormal);
    glNormal3d(SurfaceUnitNormal.x(), SurfaceUnitNormal.y(),
	       SurfaceUnitNormal.z());
    
    //Loop through the four edges of each G4Facet...
    G4bool notLastEdge;
    G4Point3D vertex[4];
    G4int edgeFlag[4], lastEdgeFlag (true);
    edgeFlag[0] = G4int (true);
    G4int edgeCount = 0;
    glEdgeFlag (GL_TRUE);
    do {
      notLastEdge = polyhedron.GetNextVertex (vertex[edgeCount], 
					      edgeFlag[edgeCount]);
      // Check to see if edge is visible or not...
      if (edgeFlag[edgeCount] != lastEdgeFlag) {
	lastEdgeFlag = edgeFlag[edgeCount];
	if (edgeFlag[edgeCount]) {
	  glEdgeFlag (GL_TRUE);
	} else {
	  glEdgeFlag (GL_FALSE);
	}
      }
      glVertex3d (vertex[edgeCount].x(), 
		  vertex[edgeCount].y(),
		  vertex[edgeCount].z());
      edgeCount++;
    } while (notLastEdge);
    while (edgeCount < 4) {  // duplicate last real vertex.
      vertex[edgeCount] = vertex[edgeCount-1];
      glVertex3d (vertex[edgeCount].x(),
		  vertex[edgeCount].y(), 
		  vertex[edgeCount].z());
      edgeCount++;
    }
    
    //do it all over again for hlr 2nd pass...
    if  (drawing_style == G4ViewParameters::hlr) {
      glEnd();
      
      //glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, clear_colour);
      glStencilFunc (GL_EQUAL, 0, 1);
      glStencilOp (GL_KEEP, GL_KEEP, GL_KEEP);
      glPolygonMode (GL_FRONT, GL_FILL);
      glColor4dv (clear_colour);
      edgeCount=0;
      
      glBegin (GL_QUADS);
      do {
	if (edgeFlag[edgeCount] != lastEdgeFlag) {
	  lastEdgeFlag = edgeFlag[edgeCount];
	  if (edgeFlag[edgeCount]) {
	    glEdgeFlag (GL_TRUE);
	  } else {
	    glEdgeFlag (GL_FALSE);
	  }
	}
	glVertex3d (vertex[edgeCount].x(), 
		    vertex[edgeCount].y(),
		    vertex[edgeCount].z());
      } while (++edgeCount < 4);
      glEnd();
      
      //and once more to reset the stencil bits...
      glStencilFunc (GL_ALWAYS, 0, 1);
      glStencilOp (GL_INVERT, GL_INVERT, GL_INVERT);
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
      edgeCount=0;
      
      glBegin(GL_QUADS);
      do {
      	if (edgeFlag[edgeCount] != lastEdgeFlag) {
      	  lastEdgeFlag = edgeFlag[edgeCount];
      	  if (edgeFlag[edgeCount]) {
      	    glEdgeFlag (GL_TRUE);
      	  } else {
      	    glEdgeFlag (GL_FALSE);
      	  }
      	}
      	glVertex3d (vertex[edgeCount].x(), 
      		    vertex[edgeCount].y(),
      		    vertex[edgeCount].z());
      } while (++edgeCount < 4);
      
    }	
    
  } while (notLastFace);  
  
  glEnd ();

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

  const G4Colour& c = GetColour (nurb);

  switch (GetDrawingStyle(nurb)) {

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

#endif
