// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLSceneHandler.cc,v 1.4 1999-12-15 14:54:08 gunter Exp $
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
  G4cerr
    << "G4OpenGLSceneHandler::AddPrimitive (const G4Text& text) not implemented yet."
    << "\n  Called with text " << text.GetText ()
    << " at " << text.GetPosition ()
    << ", size " << size
    << ", offsets " << text.GetXOffset () << ", " << text.GetYOffset ()
    << ", type " << sizeType
    << G4endl;
}

void G4OpenGLSceneHandler::AddPrimitive (const G4Circle& circle) {

  const G4Colour& c = GetColour (circle);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());
  
  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else glEnable (GL_DEPTH_TEST);
  
  glDisable (GL_LIGHTING);
  
  G4VMarker::FillStyle style = circle.GetFillStyle();
  
  switch (style) {
  case G4VMarker::noFill: 
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    break;
    
  case G4VMarker::hashed:
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glLineStipple (1, 0x0101);
    break;
    
  case G4VMarker::filled:
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;
    
  default:
    G4cerr << "Unrecognised fill style for G4Circle in G4OpenGLSceneHandler."
	 << "\nUsing G4VMarker::filled." << G4endl;
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;
    
  }
  
  G4Point3D centre = circle.GetPosition();
 
  G4bool userSpecified = (circle.GetWorldSize() || circle.GetScreenSize());
  
  const G4VMarker& def = fpViewer -> GetViewParameters().GetDefaultMarker();
  
  //Make sure we Draw circles...
  glEnable (GL_POINT_SMOOTH);

  G4double size;
  G4double scale = fpViewer -> GetViewParameters().GetGlobalMarkerScale();
  if (size = scale * // Assignment intentional.
      userSpecified ? circle.GetWorldSize() : def.GetWorldSize()) {

    // Draw in world coordinates...
    GLdouble winx, winy, winz;
    GLdouble* model_matrix = new GLdouble[16];
    GLdouble* proj_matrix = new GLdouble[16];
    GLint* v_port = new GLint[4];
    //    glPointSize (5.0); //Still have to work out correspondance between
    //                       //world and screen units.
    glGetDoublev (GL_MODELVIEW_MATRIX, model_matrix);
    glGetDoublev (GL_PROJECTION_MATRIX, proj_matrix);
    glGetIntegerv (GL_VIEWPORT, v_port);
    gluProject (size, size, size,
		model_matrix,
		proj_matrix,
		v_port,
		&winx, &winy, &winz);
    glPointSize (winx);
    glBegin (GL_POINTS);
    glVertex3d (centre.x(), centre.y(), centre.z());
    glEnd ();
    delete[] model_matrix;
    delete[] proj_matrix;
    delete[] v_port;

  }
  else {
    size = scale *
      userSpecified ? circle.GetScreenSize() : def.GetScreenSize();

    // Draw in screen coordinates...
    //    G4double* projection_matrix = new G4double[16];
    //    glGetDoublev (GL_MODELVIEW_MATRIX, projection_matrix);

    //    HepMatrix proj(4, 4, 1);
    //    for (col = 1; col < 5; col++) {
    //      for (row = 1; row < 5; row++) {
    //	proj (row, col) = projection_matrix[(col-1)*4 + (row-1)];
    //      }
    //    }

    //    G4cout << "proj is : \n" << proj << G4endl;

    //    proj.invert(ierr);

    //    G4cout << "ierr is " << ierr << G4endl;
    //    G4cout << "proj^-1 is : \n" << proj << G4endl;
    
    //    for (col = 1; col < 5; col++) {
    //      for (row = 1; row < 5; row++) {
    //	projection_matrix[(col-1)*4 + (row-1)] = proj (row, col);
    //      }
    //    }

    //    glMultMatrixd (projection_matrix);

    //    glTranslated (centre.x(), centre.y(), centre.z());
    //    glBegin (GL_LINE_LOOP);
    glPointSize (size);
    glBegin (GL_POINTS);
    //    for (G4int segment = 0; segment < num_sides; segment++) {
    //      theta = d_theta * segment;
    //      glVertex3d (sin(theta) * size,
    //		  cos(theta) * size,
    //		  0.);
    //    }
    glVertex3d (centre.x(), centre.y(), centre.z());
    glEnd ();

  }
}


void G4OpenGLSceneHandler::AddPrimitive (const G4Square& Square) {

  const G4Colour& c = GetColour (Square);
  glColor3d (c.GetRed (), c.GetGreen (), c.GetBlue ());

  if (fpViewer -> GetViewParameters ().IsMarkerNotHidden ())
    glDisable (GL_DEPTH_TEST);
  else glEnable (GL_DEPTH_TEST);

  glDisable (GL_LIGHTING);
  
  G4VMarker::FillStyle style = Square.GetFillStyle();
  
  switch (style) {
  case G4VMarker::noFill: 
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    break;
    
  case G4VMarker::hashed:
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glLineStipple (1, 0x0101);
    break;
    
  case G4VMarker::filled:
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;
    
  default:
    G4cerr << "Unrecognised fill style for G4Square in G4OpenGLSceneHandler."
	 << "\nUsing G4VMarker::filled." << G4endl;
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    break;

  }

  G4Point3D centre = Square.GetPosition();
  //  glTranslated (centre.x(), centre.y(), centre.z());

  //  G4double a = Square.GetWorldSize ();
  //  if (a <= 0.0) {
  //    a = 1000.0;
  //  }

  //  glBegin (GL_LINE_LOOP);

  //  glVertex3d ( a,  a, 0.);
  //  glVertex3d (-a,  a, 0.);
  //  glVertex3d (-a, -a, 0.);
  //  glVertex3d ( a, -a, 0.);

  //Make sure we Draw squares...
  glDisable (GL_POINT_SMOOTH);

  G4bool userSpecified = (Square.GetWorldSize() || Square.GetScreenSize());
  
  const G4VMarker& def = fpViewer -> GetViewParameters().GetDefaultMarker();
  
  G4double size;
  G4double scale = fpViewer -> GetViewParameters().GetGlobalMarkerScale();
  if (size = scale * // Assignment intentional.
      userSpecified ? Square.GetWorldSize() : def.GetWorldSize()) {

    //Draw in world coordinates...
    GLdouble winx, winy, winz;
    GLdouble* model_matrix = new GLdouble[16];
    GLdouble* proj_matrix = new GLdouble[16];
    GLint* v_port = new GLint[4];
    //    glPointSize (5.0); //Still have to work out correspondance between
    //                       //world and screen units.
    glGetDoublev (GL_MODELVIEW_MATRIX, model_matrix);
    glGetDoublev (GL_PROJECTION_MATRIX, proj_matrix);
    glGetIntegerv (GL_VIEWPORT, v_port);
    gluProject (size, size, size,
		model_matrix,
		proj_matrix,
		v_port,
		&winx, &winy, &winz);
    glPointSize (winx);
    //    glPointSize (5.0); //Still have to work out correspondance between
    //                       //world units and screen units...
    glBegin (GL_POINTS);
    glVertex3d (centre.x(), centre.y(), centre.z());
    glEnd ();
    delete[] model_matrix;
    delete[] proj_matrix;
    delete[] v_port;

  }
  else {
    size = scale *
      userSpecified ? Square.GetScreenSize() : def.GetScreenSize();

    // Draw in screen coordinates...
    glPointSize (size);
    glBegin (GL_POINTS);
    glVertex3d (centre.x(), centre.y(), centre.z());
    glEnd ();

  }
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
