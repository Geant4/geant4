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
// $Id: G4OpenGLViewer.cc,v 1.34 2007/05/24 18:27:13 allison Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL view - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4ios.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLTransform3D.hh"

#include "G4Scene.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Plane3D.hh"
#include "G4AttHolder.hh"
#include "G4AttCheck.hh"
#include <sstream>

static const char* gouraudtriangleEPS[] =
{
  "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3",
  "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd",
  "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2",
  "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse",
  "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get",
  "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get",
  "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get",
  "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll",
  "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2",
  "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4",
  "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add",
  "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1",
  "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3",
  "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll",
  "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll",
  "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array",
  "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17",
  "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index",
  "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6",
  "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index",
  "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14",
  "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3",
  "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7",
  "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add",
  "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd",
  NULL
};

G4OpenGLViewer::G4OpenGLViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
pointSize (0),
print_colour (true),
vectored_ps (true),
fOpenGLSceneHandler(scene),
background (G4Colour(0.,0.,0.)),
transparency_enabled (true),
antialiasing_enabled (false),
haloing_enabled (false),
fStartTime(-G4OPENGL_DBL_MAX),
fEndTime(G4OPENGL_DBL_MAX),
fFadeFactor(0.),
fDisplayHeadTime(false),
fDisplayHeadTimeX(-0.9),
fDisplayHeadTimeY(-0.9),
fDisplayHeadTimeSize(24.),
fDisplayHeadTimeRed(0.),
fDisplayHeadTimeGreen(1.),
fDisplayHeadTimeBlue(1.),
fDisplayLightFront(false),
fDisplayLightFrontX(0.),
fDisplayLightFrontY(0.),
fDisplayLightFrontZ(0.),
fDisplayLightFrontT(0.),
fDisplayLightFrontRed(0.),
fDisplayLightFrontGreen(1.),
fDisplayLightFrontBlue(0.)
{
  // Make changes to view parameters for OpenGL...
  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);

  //  glClearColor (0.0, 0.0, 0.0, 0.0);
  //  glClearDepth (1.0);
  //  glDisable (GL_BLEND);
  //  glDisable (GL_LINE_SMOOTH);
  //  glDisable (GL_POLYGON_SMOOTH);

  strcpy (print_string, "G4OpenGL.eps");
}

G4OpenGLViewer::~G4OpenGLViewer () {}

void G4OpenGLViewer::InitializeGLView () 
{
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glClearDepth (1.0);
  glDisable (GL_BLEND);
  glDisable (GL_LINE_SMOOTH);
  glDisable (GL_POLYGON_SMOOTH);
}  

void G4OpenGLViewer::ClearView () {
  glClearColor (background.GetRed(),
		background.GetGreen(),
		background.GetBlue(),
		1.);
  glClearDepth (1.0);
  //Below line does not compile with Mesa includes. 
  //glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
  glClear (GL_COLOR_BUFFER_BIT);
  glClear (GL_DEPTH_BUFFER_BIT);
  glClear (GL_STENCIL_BUFFER_BIT);
  glFlush ();
}

void G4OpenGLViewer::SetView () {
  
  // Calculates view representation based on extent of object being
  // viewed and (initial) viewpoint.  (Note: it can change later due
  // to user interaction via visualization system's GUI.)
  
  // Lighting.
  GLfloat lightPosition [4];
  lightPosition [0] = fVP.GetActualLightpointDirection().x();
  lightPosition [1] = fVP.GetActualLightpointDirection().y();
  lightPosition [2] = fVP.GetActualLightpointDirection().z();
  lightPosition [3] = 0.;
  // Light position is "true" light direction, so must come after gluLookAt.
  GLfloat ambient [] = { 0.2, 0.2, 0.2, 1.};
  GLfloat diffuse [] = { 0.8, 0.8, 0.8, 1.};
  glEnable (GL_LIGHT0);
  glLightfv (GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, diffuse);
  
  // Get radius of scene, etc.
  // Note that this procedure properly takes into account zoom, dolly and pan.
  const G4Point3D targetPoint
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Point3D cameraPosition =
    targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
  const GLdouble pnear   = fVP.GetNearDistance (cameraDistance, radius);
  const GLdouble pfar    = fVP.GetFarDistance  (cameraDistance, pnear, radius);
  const GLdouble right  = fVP.GetFrontHalfHeight (pnear, radius);
  const GLdouble left   = -right;
  const GLdouble bottom = left;
  const GLdouble top    = right;
  
  glMatrixMode (GL_PROJECTION); // set up Frustum.
  glLoadIdentity();

  const G4Vector3D scale = fVP.GetScaleFactor();
  glScaled(scale.x(),scale.y(),scale.z());
  
  if (fVP.GetFieldHalfAngle() == 0.) {
    glOrtho (left, right, bottom, top, pnear, pfar);
  }
  else {
    glFrustum (left, right, bottom, top, pnear, pfar);
  }
  
  glMatrixMode (GL_MODELVIEW); // apply further transformations to scene.
  glLoadIdentity();
  
  const G4Normal3D& upVector = fVP.GetUpVector ();  
  G4Point3D gltarget;
  if (cameraDistance > 1.e-6 * radius) {
    gltarget = targetPoint;
  }
  else {
    gltarget = targetPoint - radius * fVP.GetViewpointDirection().unit();
  }

  const G4Point3D& pCamera = cameraPosition;  // An alias for brevity.
  gluLookAt (pCamera.x(),  pCamera.y(),  pCamera.z(),       // Viewpoint.
	     gltarget.x(), gltarget.y(), gltarget.z(),      // Target point.
	     upVector.x(), upVector.y(), upVector.z());     // Up vector.
  
  // Light position is "true" light direction, so must come after gluLookAt.
  glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

  // OpenGL no longer seems to reconstruct clipped edges, so, when the
  // BooleanProcessor is up to it, abandon this and use generic
  // clipping in G4OpenGLSceneHandler::CreateSectionPolyhedron.  Also,
  // force kernel visit on change of clipping plane in
  // G4OpenGLStoredViewer::CompareForKernelVisit.
  if (fVP.IsSection () ) {  // pair of back to back clip planes.
    const G4Plane3D& s = fVP.GetSectionPlane ();
    double sArray[4];
    sArray[0] = s.a();
    sArray[1] = s.b();
    sArray[2] = s.c();
    sArray[3] = s.d() + radius * 1.e-05;
    glClipPlane (GL_CLIP_PLANE0, sArray);
    glEnable (GL_CLIP_PLANE0);
    sArray[0] = -s.a();
    sArray[1] = -s.b();
    sArray[2] = -s.c();
    sArray[3] = -s.d() + radius * 1.e-05;
    glClipPlane (GL_CLIP_PLANE1, sArray);
    glEnable (GL_CLIP_PLANE1);
  } else {
    glDisable (GL_CLIP_PLANE0);
    glDisable (GL_CLIP_PLANE1);
  }

  const G4Planes& cutaways = fVP.GetCutawayPlanes();
  size_t nPlanes = cutaways.size();
  if (fVP.IsCutaway() &&
      fVP.GetCutawayMode() == G4ViewParameters::cutawayIntersection &&
      nPlanes > 0) {
    double a[4];
    a[0] = cutaways[0].a();
    a[1] = cutaways[0].b();
    a[2] = cutaways[0].c();
    a[3] = cutaways[0].d();
    glClipPlane (GL_CLIP_PLANE2, a);
    glEnable (GL_CLIP_PLANE2);
    if (nPlanes > 1) {
      a[0] = cutaways[1].a();
      a[1] = cutaways[1].b();
      a[2] = cutaways[1].c();
      a[3] = cutaways[1].d();
      glClipPlane (GL_CLIP_PLANE3, a);
      glEnable (GL_CLIP_PLANE3);
    }
    if (nPlanes > 2) {
      a[0] = cutaways[2].a();
      a[1] = cutaways[2].b();
      a[2] = cutaways[2].c();
      a[3] = cutaways[2].d();
      glClipPlane (GL_CLIP_PLANE4, a);
      glEnable (GL_CLIP_PLANE4);
    }
  } else {
    glDisable (GL_CLIP_PLANE2);
    glDisable (GL_CLIP_PLANE3);
    glDisable (GL_CLIP_PLANE4);
  }

  // Background.
  background = fVP.GetBackgroundColour ();

}

void G4OpenGLViewer::HaloingFirstPass () {
  
  //To perform haloing, first Draw all information to the depth buffer
  //alone, using a chunky line width, and then Draw all info again, to
  //the colour buffer, setting a thinner line width an the depth testing 
  //function to less than or equal, so if two lines cross, the one 
  //passing behind the other will not pass the depth test, and so not
  //get rendered either side of the infront line for a short distance.

  //First, disable writing to the colo(u)r buffer...
  glColorMask (GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

  //Now enable writing to the depth buffer...
  glDepthMask (GL_TRUE);
  glDepthFunc (GL_LESS);
  glClearDepth (1.0);

  //Finally, set the line width to something wide...
  glLineWidth (3.0);

}

void G4OpenGLViewer::HaloingSecondPass () {

  //And finally, turn the colour buffer back on with a sesible line width...
  glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDepthFunc (GL_LEQUAL);
  glLineWidth (1.0);

}

void G4OpenGLViewer::Pick(GLdouble x, GLdouble y)
{
  //G4cout << "X: " << x << ", Y: " << y << G4endl;
  const G4int BUFSIZE = 512;
  GLuint selectBuffer[BUFSIZE];
  glSelectBuffer(BUFSIZE, selectBuffer);
  glRenderMode(GL_SELECT);
  glInitNames();
  glPushName(0);
  glMatrixMode(GL_PROJECTION);
  G4double currentProjectionMatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX, currentProjectionMatrix);
  glPushMatrix();
  glLoadIdentity();
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  // Define 5x5 pixel pick area
  gluPickMatrix(x, viewport[3] - y, 5., 5., viewport);
  glMultMatrixd(currentProjectionMatrix);
  glMatrixMode(GL_MODELVIEW);
  DrawView();
  GLint hits = glRenderMode(GL_RENDER);
  if (hits < 0)
    G4cout << "Too many hits.  Zoom in to reduce overlaps." << G4cout;
  else if (hits > 0) {
    //G4cout << hits << " hit(s)" << G4endl;
    GLuint* p = selectBuffer;
    for (GLint i = 0; i < hits; ++i) {
      GLuint nnames = *p++;
      *p++; //OR GLuint zmin = *p++;
      *p++; //OR GLuint zmax = *p++;
      //G4cout << "Hit " << i << ": " << nnames << " names"
      //     << "\nzmin: " << zmin << ", zmax: " << zmax << G4endl;
      for (GLuint j = 0; j < nnames; ++j) {
	GLuint name = *p++;
	//G4cout << "Name " << j << ": PickName: " << name << G4endl;
	std::map<GLuint, G4AttHolder*>::iterator iter =
	  fOpenGLSceneHandler.fPickMap.find(name);
	if (iter != fOpenGLSceneHandler.fPickMap.end()) {
	  G4AttHolder* attHolder = iter->second;
	  if(attHolder && attHolder->GetAttDefs().size()) {
	    for (size_t i = 0; i < attHolder->GetAttDefs().size(); ++i) {
	      G4cout << G4AttCheck(attHolder->GetAttValues()[i],
				   attHolder->GetAttDefs()[i]);
	    }
	  }
	}
      }
      G4cout << G4endl;
    }
  }
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

void G4OpenGLViewer::print() {

  // Print vectored PostScript
  
  G4int size = 5000000;

  GLfloat* feedback_buffer;
  GLint returned;
  FILE* file;
  
  feedback_buffer = new GLfloat[size];
  glFeedbackBuffer (size, GL_3D_COLOR, feedback_buffer);
  glRenderMode (GL_FEEDBACK);
  
  DrawView();
  returned = glRenderMode (GL_RENDER);
  
  if (print_string) {
    file = fopen (print_string, "w");
    if (file) {
      spewWireframeEPS (file, returned, feedback_buffer, "rendereps");
    } else {
      printf("Could not open %s\n", print_string);
    }
  } else {
    printBuffer (returned, feedback_buffer);
  }

  delete[] feedback_buffer;
}

void G4OpenGLViewer::print3DcolorVertex(GLint size, GLint * count, GLfloat * buffer)
{
  G4int i;

  printf("  ");
  for (i = 0; i < 7; i++) {
    printf("%4.2f ", buffer[size - (*count)]);
    *count = *count - 1;
  }
  printf("\n");
}

void G4OpenGLViewer::spewWireframeEPS (FILE* file, GLint size, GLfloat* buffer, const char* cr) {

  GLfloat EPS_GOURAUD_THRESHOLD=0.1;

  GLfloat clearColor[4], viewport[4];
  GLfloat lineWidth;
  G4int i;

  glGetFloatv (GL_VIEWPORT, viewport);
  glGetFloatv (GL_COLOR_CLEAR_VALUE, clearColor);
  glGetFloatv (GL_LINE_WIDTH, &lineWidth);
  glGetFloatv (GL_POINT_SIZE, &pointSize);

  fputs ("%!PS-Adobe-2.0 EPSF-2.0\n", file);
  fprintf (file, "%%%%Creator: %s (using OpenGL feedback)\n", cr);
  fprintf (file, "%%%%BoundingBox: %g %g %g %g\n", viewport[0], viewport[1], viewport[2], viewport[3]);
  fputs ("%%EndComments\n", file);
  fputs ("\n", file);
  fputs ("gsave\n", file);
  fputs ("\n", file);

  fputs ("% the gouraudtriangle PostScript fragment below is free\n", file);
  fputs ("% written by Frederic Delhoume (delhoume@ilog.fr)\n", file);
  fprintf (file, "/threshold %g def\n", EPS_GOURAUD_THRESHOLD);
  for (i=0; gouraudtriangleEPS[i]; i++) {
    fprintf (file, "%s\n", gouraudtriangleEPS[i]);
  }

  fprintf(file, "\n%g setlinewidth\n", lineWidth);
  
  fprintf (file, "%g %g %g setrgbcolor\n", clearColor[0], clearColor[1], clearColor[2]);
  fprintf (file, "%g %g %g %g rectfill\n\n", viewport[0], viewport[1], viewport[2], viewport[3]);

  spewSortedFeedback (file, size, buffer);

  fputs ("grestore\n\n", file);
  fputs ("showpage\n", file);

  fclose(file);
}

void G4OpenGLViewer::printBuffer (GLint size, GLfloat* buffer) {

  GLint count;
  G4int token, nvertices;

  count=size;
  while(count) {
    token=G4int (buffer[size-count]);
    count--;
    switch (token) {

    case GL_PASS_THROUGH_TOKEN:
      printf ("GL_PASS_THROUGH_TOKEN\n");
      printf ("  %4.2f\n", buffer[size-count]);
      count--;
      break;

    case GL_POINT_TOKEN:
      printf ("GL_POINT_TOKEN\n");
      print3DcolorVertex (size, &count, buffer);
      break;

    case GL_LINE_TOKEN:
      printf ("GL_LINE_TOKEN\n");
      print3DcolorVertex (size, &count, buffer);
      print3DcolorVertex (size, &count, buffer);
      break;
      
    case GL_LINE_RESET_TOKEN:
      printf ("GL_LINE_RESET_TOKEN\n");
      print3DcolorVertex (size, &count, buffer);
      print3DcolorVertex (size, &count, buffer);
      break;

    case GL_POLYGON_TOKEN:
      printf ("GL_POLYGON_TOKEN\n");
      nvertices=G4int (buffer[size-count]);
      count--;
      for (; nvertices>0; nvertices--) {
	print3DcolorVertex (size, &count, buffer);
      }
    }
  }
}

G4float* G4OpenGLViewer::spewPrimitiveEPS (FILE* file, GLfloat* loc) {
  
  G4int token;
  G4int nvertices, i;
  GLfloat red, green, blue, intensity;
  G4int smooth;
  GLfloat dx, dy, dr, dg, db, absR, absG, absB, colormax;
  G4int steps;
  Feedback3Dcolor *vertex;
  GLfloat xstep(0.), ystep(0.), rstep(0.), gstep(0.), bstep(0.);
  GLfloat xnext(0.), ynext(0.), rnext(0.), gnext(0.), bnext(0.), distance(0.);

  token=G4int (*loc);
  loc++;
  switch (token) {
  case GL_LINE_RESET_TOKEN:
  case GL_LINE_TOKEN:
    vertex=(Feedback3Dcolor*)loc;
    dr=vertex[1].red - vertex[0].red;
    dg=vertex[1].green - vertex[0].green;
    db=vertex[1].blue - vertex[0].blue;

    if (!print_colour) {
      dr+=(dg+db);
      dr/=3.0;
      dg=dr;
      db=dr;
    }

    if (dr!=0 || dg!=0 || db!=0) {
      dx=vertex[1].x - vertex[0].x;
      dy=vertex[1].y - vertex[0].y;
      distance=std::sqrt(dx*dx + dy*dy);

      absR=std::fabs(dr);
      absG=std::fabs(dg);
      absB=std::fabs(db);

      #define Max(a, b) (((a)>(b))?(a):(b))

      #define EPS_SMOOTH_LINE_FACTOR 0.06

      colormax=Max(absR, Max(absG, absB));
      steps=Max(1, G4int (colormax*distance*EPS_SMOOTH_LINE_FACTOR));
      
      xstep=dx/steps;
      ystep=dy/steps;

      rstep=dr/steps;
      gstep=dg/steps;
      bstep=db/steps;

      xnext=vertex[0].x;
      ynext=vertex[0].y;
      rnext=vertex[0].red;
      gnext=vertex[0].green;
      bnext=vertex[0].blue;

      if (!print_colour) {
	rnext+=(gnext+bnext);
	rnext/=3.0;
	gnext=rnext;
	bnext=rnext;
      }

      xnext -= xstep/2.0;
      ynext -= ystep/2.0;
      rnext -= rstep/2.0;
      gnext -= gstep/2.0;
      bnext -= bstep/2.0;
    } else {
      steps=0;
    }
    if (print_colour) {
      fprintf (file, "%g %g %g setrgbcolor\n",
	       vertex[0].red, vertex[0].green, vertex[0].blue);
    } else {
      intensity = (vertex[0].red + vertex[0].green + vertex[0].blue) / 3.0;
      fprintf (file, "%g %g %g setrgbcolor\n",
	       intensity, intensity, intensity);
    }      
    fprintf (file, "%g %g moveto\n", vertex[0].x, vertex[0].y);

    for (i=0; i<steps; i++) {

      xnext += xstep;
      ynext += ystep;
      rnext += rstep;
      gnext += gstep;
      bnext += bstep;

      fprintf (file, "%g %g lineto stroke\n", xnext, ynext);
      fprintf (file, "%g %g %g setrgbcolor\n", rnext, gnext, bnext);
      fprintf (file, "%g %g moveto\n", xnext, ynext);
    }
    fprintf (file, "%g %g lineto stroke\n", vertex[1].x, vertex[1].y);

    loc += 14;
    break;

  case GL_POLYGON_TOKEN:
    nvertices = G4int (*loc);
    loc++;
    vertex=(Feedback3Dcolor*)loc;
    if (nvertices>0) {
      red=vertex[0].red;
      green=vertex[0].green;
      blue=vertex[0].blue;
      smooth=0;
      
      if (!print_colour) {
	red+=(green+blue);
	red/=3.0;
	green=red;
	blue=red;
      }
      
      if (print_colour) {
	for (i=1; i<nvertices; i++) {
	  if (red!=vertex[i].red || green!=vertex[i].green || blue!=vertex[i].blue) {
	    smooth=1;
	    break;
	  }
	}
      } else {
	for (i=1; i<nvertices; i++) {
	  intensity = vertex[i].red + vertex[i].green + vertex[i].blue;
	  intensity/=3.0;
	  if (red!=intensity) {
	    smooth=1;
	    break;
	  }
	}
      }

      if (smooth) {
	G4int triOffset;
	for (i=0; i<nvertices-2; i++) {
	  triOffset = i*7;
	  fprintf (file, "[%g %g %g %g %g %g]",
		   vertex[0].x, vertex[i+1].x, vertex[i+2].x,
		   vertex[0].y, vertex[i+1].y, vertex[i+2].y);
	  if (print_colour) {
	    fprintf (file, " [%g %g %g] [%g %g %g] [%g %g %g] gouraudtriangle\n",
		     vertex[0].red, vertex[0].green, vertex[0].blue,
		     vertex[i+1].red, vertex[i+1].green, vertex[i+1].blue,
		     vertex[i+2].red, vertex[i+2].green, vertex[i+2].blue);
	  } else {

	    intensity = vertex[0].red + vertex[0].green + vertex[0].blue;
	    intensity/=3.0;
	    fprintf (file, " [%g %g %g]", intensity, intensity, intensity);

	    intensity = vertex[1].red + vertex[1].green + vertex[1].blue;
	    intensity/=3.0;
	    fprintf (file, " [%g %g %g]", intensity, intensity, intensity);

	    intensity = vertex[2].red + vertex[2].green + vertex[2].blue;
	    intensity/=3.0;
	    fprintf (file, " [%g %g %g] gouraudtriangle\n", intensity, intensity, intensity);
	  }
	}
      } else {
	fprintf (file, "newpath\n");
	fprintf (file, "%g %g %g setrgbcolor\n", red, green, blue);
	fprintf (file, "%g %g moveto\n", vertex[0].x, vertex[0].y);
	for (i=1; i<nvertices; i++) {
	  fprintf (file, "%g %g lineto\n", vertex[i].x, vertex[i].y);
	}
	fprintf (file, "closepath fill\n\n");
      }
    }
    loc += nvertices*7;
    break;

  case GL_POINT_TOKEN:
    vertex=(Feedback3Dcolor*)loc;
    if (print_colour) {
      fprintf (file, "%g %g %g setrgbcolor\n", vertex[0].red, vertex[0].green, vertex[0].blue);
    } else {
      intensity = vertex[0].red + vertex[0].green + vertex[0].blue;
      intensity/=3.0;
      fprintf (file, "%g %g %g setrgbcolor\n", intensity, intensity, intensity);
    }      
    fprintf(file, "%g %g %g 0 360 arc fill\n\n", vertex[0].x, vertex[0].y, pointSize / 2.0);
    loc += 7;           /* Each vertex element in the feedback
                           buffer is 7 GLfloats. */
    break;
  default:
    /* XXX Left as an excersie to the reader. */
    static G4bool spewPrimitiveEPSWarned = false;
    if (!spewPrimitiveEPSWarned) {
      std::ostringstream oss;
      oss <<
	"Incomplete implementation.  Unexpected token (" << token << ")."
	"\n  (Seems to be caused by text.)";
      G4Exception("G4OpenGLViewer::spewPrimitiveEPS",
		  "Unexpected token",
		  JustWarning,
		  oss.str().c_str());
      spewPrimitiveEPSWarned = true;
    }
  }
  return loc;
}

typedef struct G4OpenGLViewerDepthIndex {
  GLfloat *ptr;
  GLfloat depth;
} DepthIndex;

extern "C" {
  int G4OpenGLViewercompare(const void *a, const void *b)
  {
    const DepthIndex *p1 = (DepthIndex *) a;
    const DepthIndex *p2 = (DepthIndex *) b;
    GLfloat diff = p2->depth - p1->depth;
    
    if (diff > 0.0) {
      return 1;
    } else if (diff < 0.0) {
      return -1;
    } else {
      return 0;
    }
  }
}

void G4OpenGLViewer::spewSortedFeedback(FILE * file, GLint size, GLfloat * buffer)
{
  int token;
  GLfloat *loc, *end;
  Feedback3Dcolor *vertex;
  GLfloat depthSum;
  int nprimitives, item;
  DepthIndex *prims;
  int nvertices, i;

  end = buffer + size;

  /* Count how many primitives there are. */
  nprimitives = 0;
  loc = buffer;
  while (loc < end) {
    token = int (*loc);
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      loc += 14;
      nprimitives++;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = int (*loc);
      loc++;
      loc += (7 * nvertices);
      nprimitives++;
      break;
    case GL_POINT_TOKEN:
      loc += 7;
      nprimitives++;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      static G4bool spewSortedFeedbackWarned = false;
      if (!spewSortedFeedbackWarned) {
	std::ostringstream oss;
	oss <<
	  "Incomplete implementation.  Unexpected token (" << token << ")."
	  "\n  (Seems to be caused by text.)";
	G4Exception("G4OpenGLViewer::spewSortedFeedback",
		    "Unexpected token",
		    JustWarning,
		    oss.str().c_str());
	spewSortedFeedbackWarned = true;
      }
      nprimitives++;
    }
  }

  /* Allocate an array of pointers that will point back at
     primitives in the feedback buffer.  There will be one
     entry per primitive.  This array is also where we keep the
     primitive's average depth.  There is one entry per
     primitive  in the feedback buffer. */
  prims = (DepthIndex *) malloc(sizeof(DepthIndex) * nprimitives);

  item = 0;
  loc = buffer;
  while (loc < end) {
    prims[item].ptr = loc;  /* Save this primitive's location. */
    token = int (*loc);
    loc++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z + vertex[1].z;
      prims[item].depth = depthSum / 2.0;
      loc += 14;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = int (*loc);
      loc++;
      vertex = (Feedback3Dcolor *) loc;
      depthSum = vertex[0].z;
      for (i = 1; i < nvertices; i++) {
        depthSum += vertex[i].z;
      }
      prims[item].depth = depthSum / nvertices;
      loc += (7 * nvertices);
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor *) loc;
      prims[item].depth = vertex[0].z;
      loc += 7;
      break;
    default:
      /* XXX Left as an excersie to the reader. */
      assert(1);
    }
    item++;
  }
  assert(item == nprimitives);

  /* Sort the primitives back to front. */
  qsort(prims, nprimitives, sizeof(DepthIndex), G4OpenGLViewercompare);

  /* Understand that sorting by a primitives average depth
     doesn't allow us to disambiguate some cases like self
     intersecting polygons.  Handling these cases would require
     breaking up the primitives.  That's too involved for this
     example.  Sorting by depth is good enough for lots of
     applications. */

  /* Emit the Encapsulated PostScript for the primitives in
     back to front order. */
  for (item = 0; item < nprimitives; item++) {
    (void) spewPrimitiveEPS(file, prims[item].ptr);
  }

  free(prims);
}

#endif
