// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXViewer.cc,v 1.9 2001-02-03 18:39:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// G4OpenGLXViewer : Class to provide XWindows specific
//                 functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLXViewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>
#include <string.h>

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

#include <X11/Xatom.h>
#include <X11/Xutil.h>

#include <stdlib.h>
#include <string.h>

#define NewString(str) \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)

#define USE_DEFAULT_COLORMAP 1
#define USE_STANDARD_COLORMAP 0

XVisualInfo*  G4OpenGLXViewer::vi_single_buffer = 0;
XVisualInfo*  G4OpenGLXViewer::vi_double_buffer = 0;

extern "C"
{
  static Bool WaitForNotify (Display*, XEvent* e, char* arg) {
    return (e->type == MapNotify) && (e->xmap.window == (Window) arg);
  }
}

void G4OpenGLXViewer::SetView () {
  glXMakeCurrent (dpy, win, cx);
  G4OpenGLViewer::SetView ();  
}

void G4OpenGLXViewer::ShowView () {
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propogated before progressing.
  glFlush ();
}

void G4OpenGLXViewer::FinishView () {
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propogated before progressing.
  if (doublebuffer == true) {
    glXSwapBuffers (dpy, win);  
  }
  else glFlush ();
}

void G4OpenGLXViewer::GetXConnection () {
// get a connection.
  dpy = XOpenDisplay (0);
  if (!dpy) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer couldn't open display." << G4endl;
    return;
  }

// make sure OpenGL is supported and installed properly.
  if (!glXQueryExtension (dpy, &errorBase, &eventBase)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer X Server has no GLX extension." 
	 << G4endl;
    return;
  }
}

void G4OpenGLXViewer::CreateGLXContext (XVisualInfo* v) {

  vi = v;
// get window's attributes
  if (!XGetWindowAttributes(dpy, XRootWindow (dpy, vi -> screen), &xwa)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer couldn't return window attributes"
	 << G4endl;
    return;
  }
  
// create a GLX context
  cx = glXCreateContext (dpy, vi, 0, true);
  if (!cx) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer couldn't create context."
	 << G4endl;
    return;
  }

// New stab at getting a colormap

  Status status;
  XStandardColormap *standardCmaps = XAllocStandardColormap ();
  int i, numCmaps;

  status = XmuLookupStandardColormap (dpy, 
				      vi -> screen, 
				      vi -> visualid, 
				      vi -> depth, 
				      XA_RGB_DEFAULT_MAP, 
				      False, 
				      True);
  
  if (status == 1) {
    status = XGetRGBColormaps (dpy, 
			       XRootWindow (dpy, vi -> screen), 
			       &standardCmaps, 
			       &numCmaps, 
			       XA_RGB_DEFAULT_MAP);
    if (status == 1)
      for (i = 0; i < numCmaps; i++)
	if (standardCmaps[i].visualid == vi -> visualid) {
	  cmap = standardCmaps[i].colormap;
	  XFree (standardCmaps);
	}
    G4cout << "Got standard cmap" << G4endl;
  } else {
    cmap = XCreateColormap (dpy, 
			    XRootWindow(dpy, vi -> screen), 
			    vi -> visual, 
			    AllocNone);
    G4cout << "Created own cmap" << G4endl;
  }

  if (!cmap) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer failed to allocate a Colormap."
	 << G4endl;
    return;
  }

}
  
void G4OpenGLXViewer::CreateMainWindow () {
  
// create a window
  swa.colormap = cmap;
  swa.border_pixel = 0;
  swa.event_mask = ExposureMask | ButtonPressMask | StructureNotifyMask;
  swa.backing_store = WhenMapped;

  //  G4int                             WinSize_x;
  //  G4int                             WinSize_y;
  if (xwa.width > xwa.height) {
    WinSize_x = (xwa.height)/2;
    WinSize_y = (xwa.height)/2;
  }
  else {
    WinSize_x = (xwa.width)/2;
    WinSize_y = (xwa.width)/2;
  }
  if (WinSize_x < fVP.GetWindowSizeHintX ())
    WinSize_x = fVP.GetWindowSizeHintX ();
  if (WinSize_y < fVP.GetWindowSizeHintY ())
    WinSize_y = fVP.GetWindowSizeHintY ();
  x_origin = xwa.x;
  y_origin = xwa.y;
  G4cout << "Window name: " << fName << G4endl;
  strncpy (charViewName, fName, 100);
  char *window_name = charViewName;
  char *icon_name = charViewName;
  char tmpatom[] = "XA_WM_NORMAL_HINTS"; 
  size_hints = XAllocSizeHints();
  wm_hints = XAllocWMHints();
  class_hints = XAllocClassHint();

  size_hints -> flags = PPosition | PSize | PMinSize;
  size_hints -> min_width = 300;
  size_hints -> min_height = 200;

  XStringListToTextProperty (&window_name, 1, &windowName);
  XStringListToTextProperty (&icon_name, 1, &iconName);

  wm_hints -> initial_state = NormalState;
  wm_hints -> input = True;
  wm_hints -> icon_pixmap = icon_pixmap;
  wm_hints -> flags = StateHint | IconPixmapHint | InputHint;

  class_hints -> res_name  = NewString("G4OpenGL");
  class_hints -> res_class = NewString("G4OpenGL");

  win = XCreateWindow (dpy, XRootWindow (dpy, vi -> screen), x_origin, 
                       y_origin, WinSize_x, WinSize_y, 0, vi -> depth,
                       InputOutput, vi -> visual,  
                       CWBorderPixel | CWColormap | 
                       CWEventMask | CWBackingStore,
                       &swa);
  
  XSetWMProperties (dpy, win, &windowName, &iconName, 0, 0, 
                    size_hints, wm_hints, class_hints);
  
// request X to Draw window on screen.
  XMapWindow (dpy, win);

// Wait for window to appear (wait for an "expose" event).
  XIfEvent (dpy, &event, WaitForNotify, (char*) win);

// connect the context to a window
  glXMakeCurrent (dpy, win, cx);
}

G4OpenGLXViewer::G4OpenGLXViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene), //ajw
vi_immediate (0),
vi_stored (0),
print_colour (true),
vectored_ps (true)
{

  strcpy (print_string, "G4OpenGL.eps");

  GetXConnection ();
  if (fViewId < 0) return;
  
  // Try for a visual suitable for OpenGLImmediate..
  // first try for a single buffered RGB window
  if (!vi_single_buffer) {
    vi_single_buffer =
      glXChooseVisual (dpy, XDefaultScreen (dpy), snglBuf_RGBA);
  }
  if (!vi_single_buffer) {
    G4cout <<
    "G4OpenGLXViewer::G4OpenGLXViewer: unable to get a single buffer visual."
	   << G4endl;
  }

  if (!vi_double_buffer) {
    vi_double_buffer =
      glXChooseVisual (dpy, XDefaultScreen (dpy), dblBuf_RGBA);
  }
  if (!vi_double_buffer) {
    G4cout <<
    "G4OpenGLXViewer::G4OpenGLXViewer: unable to get a double buffer visual."
	   << G4endl;
  }

  if (vi_single_buffer) {
    vi_immediate = vi_single_buffer;
    attributeList = snglBuf_RGBA;
    doublebuffer = false;
  }
  
  if (!vi_immediate){
    // next try for a double buffered RGB, but Draw to top buffer
    if (vi_double_buffer) {
      vi_immediate = vi_double_buffer;
      attributeList = dblBuf_RGBA;
      doublebuffer = true;
    }
  }

  // Now try for a visual suitable for OpenGLStored...
  // Try for a double buffered RGB window
  if (vi_double_buffer) {
    vi_stored = vi_double_buffer;
    attributeList = dblBuf_RGBA;
    doublebuffer = true;
  }

  if (!vi_immediate || !vi_stored) {
    G4cout <<
    "G4OpenGLXViewer::G4OpenGLXViewer: unable to get required visuals."
	   << G4endl;
    fViewId = -1;  // This flags an error.
  }

  //  glClearColor (0., 0., 0., 0.);
  //  glClearDepth (1.);
  
}

G4OpenGLXViewer::~G4OpenGLXViewer () {
  if (fViewId >= 0) {
    //Close a window from here
    glXDestroyContext (dpy, cx);
    glXMakeCurrent (dpy, None, NULL);
    if (win) XDestroyWindow (dpy, win); // ...if already deleted in
    // sub-class G4OpenGLXmViewer.
    XFlush (dpy);
  }
}

void G4OpenGLXViewer::print() {
  
  //cout << "print_col_callback requested with file name: " << print_string << G4endl;
  
  if (vectored_ps) {
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
    //  free (feedback_buffer);
    delete[] feedback_buffer;

  } else {

    XVisualInfo* pvi;
    G4bool new_db;
    G4bool last_db;
    GLXContext pcx = create_GL_print_context(pvi, new_db);
    GLXContext tmp_cx;
    tmp_cx = cx;
    cx=pcx;
    last_db=doublebuffer;
    doublebuffer=new_db;
    
    Pixmap pmap = XCreatePixmap (dpy,
				 XRootWindow (dpy, pvi->screen),
				 WinSize_x, WinSize_y,
				 pvi->depth);
    
    GLXPixmap glxpmap = glXCreateGLXPixmap (dpy, 
					    pvi,
					    pmap);
    
    GLXDrawable tmp_win;
    tmp_win=win;
    win=glxpmap;
    
    glXMakeCurrent (dpy,
		    glxpmap,
		    cx);
    
    glViewport (0, 0, WinSize_x, WinSize_y);
    
    ClearView ();
    DrawView ();
    
    generateEPS (print_string,
		 print_colour,
		 WinSize_x, WinSize_y);
    
    win=tmp_win;
    cx=tmp_cx;
    doublebuffer=last_db;
    
    glXMakeCurrent (dpy,
		    glxpmap,
		    cx);
    
  }

}

void G4OpenGLXViewer::print3DcolorVertex(GLint size, GLint * count, GLfloat * buffer)
{
  G4int i;

  printf("  ");
  for (i = 0; i < 7; i++) {
    printf("%4.2f ", buffer[size - (*count)]);
    *count = *count - 1;
  }
  printf("\n");
}

void G4OpenGLXViewer::spewWireframeEPS (FILE* file, GLint size, GLfloat* buffer, const char* cr) {

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

void G4OpenGLXViewer::printBuffer (GLint size, GLfloat* buffer) {

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

G4float* G4OpenGLXViewer::spewPrimitiveEPS (FILE* file, GLfloat* loc) {
  
  G4int token;
  G4int nvertices, i;
  GLfloat red, green, blue, intensity;
  G4int smooth;
  GLfloat dx, dy, dr, dg, db, absR, absG, absB, colormax;
  G4int steps;
  Feedback3Dcolor *vertex;
  GLfloat xstep, ystep, rstep, gstep, bstep;
  GLfloat xnext, ynext, rnext, gnext, bnext, distance;

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
      distance=sqrt(dx*dx + dy*dy);

      absR=fabs(dr);
      absG=fabs(dg);
      absB=fabs(db);

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
    printf("Incomplete implementation.  Unexpected token (%d).\n", token);
    exit(1);
  }
  return loc;
}

typedef struct _DepthIndex {
  GLfloat *ptr;
  GLfloat depth;
} DepthIndex;

extern "C"
{
  static int
  compare(const void *a, const void *b)
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

void G4OpenGLXViewer::spewSortedFeedback(FILE * file, GLint size, GLfloat * buffer)
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
      printf("Incomplete implementation.  Unexpected token (%d).\n",
        token);
      exit(1);
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
  qsort(prims, nprimitives, sizeof(DepthIndex), compare);

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

GLXContext G4OpenGLXViewer::create_GL_print_context(XVisualInfo*& pvi, G4bool& db) {
  
  pvi = glXChooseVisual (dpy,
			 XDefaultScreen (dpy),
			 snglBuf_RGBA);

  if (pvi != NULL) {
    db=false;
  } else {
    pvi = glXChooseVisual (dpy,
			   XDefaultScreen (dpy),
			   dblBuf_RGBA);
    if (NULL) {
      db=true;
    }
  }

  return glXCreateContext (dpy,
			   pvi,
			   NULL,
			   False);
}

int G4OpenGLXViewer::generateEPS (char* filnam,
				int inColour,
				unsigned int width,
				unsigned int height) {

  FILE* fp;
  GLubyte* pixels;
  GLubyte* curpix;
  int components, pos, i;

  pixels = grabPixels (inColour, width, height);

  if (pixels == NULL)
    return 1;
  if (inColour) {
    components = 3;
  } else {
    components = 1;
  }
  
  fp = fopen (filnam, "w");
  if (fp == NULL) {
    return 2;
  }
  
  fprintf (fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf (fp, "%%%%Title: %s\n", filnam);
  fprintf (fp, "%%%%Creator: OpenGL pixmap render output\n");
  fprintf (fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
  fprintf (fp, "%%%%EndComments\n");
  fprintf (fp, "gsave\n");
  fprintf (fp, "/bwproc {\n");
  fprintf (fp, "    rgbproc\n");
  fprintf (fp, "    dup length 3 idiv string 0 3 0 \n");
  fprintf (fp, "    5 -1 roll {\n");
  fprintf (fp, "    add 2 1 roll 1 sub dup 0 eq\n");
  fprintf (fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
  fprintf (fp, "       3 1 roll 5 -1 roll } put 1 add 3 0 \n");
  fprintf (fp, "    { 2 1 roll } ifelse\n");
  fprintf (fp, "    }forall\n");
  fprintf (fp, "    pop pop pop\n");
  fprintf (fp, "} def\n");
  fprintf (fp, "systemdict /colorimage known not {\n");
  fprintf (fp, "   /colorimage {\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       /rgbproc exch def\n");
  fprintf (fp, "       { bwproc } image\n");
  fprintf (fp, "   }  def\n");
  fprintf (fp, "} if\n");
  fprintf (fp, "/picstr %d string def\n", width * components);
  fprintf (fp, "%d %d scale\n", width, height);
  fprintf (fp, "%d %d %d\n", width, height, 8);
  fprintf (fp, "[%d 0 0 %d 0 0]\n", width, height);
  fprintf (fp, "{currentfile picstr readhexstring pop}\n");
  fprintf (fp, "false %d\n", components);
  fprintf (fp, "colorimage\n");
  
  curpix = (GLubyte*) pixels;
  pos = 0;
  for (i = width*height*components; i>0; i--) {
    fprintf (fp, "%02hx ", *(curpix++));
    if (++pos >= 32) {
      fprintf (fp, "\n");
      pos = 0; 
    }
  }
  if (pos)
    fprintf (fp, "\n");

  fprintf (fp, "grestore\n");
  fprintf (fp, "showpage\n");
  delete pixels;
  fclose (fp);
  return 0;
}

GLubyte* G4OpenGLXViewer::grabPixels (int inColor, unsigned int width, unsigned int height) {
  
  GLubyte* buffer;
  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  GLenum format;
  int size;

  if (inColor) {
    format = GL_RGB;
    size = width*height*3;
  } else {
    format = GL_LUMINANCE;
    size = width*height*1;
  }

  buffer = new GLubyte[size];
  if (buffer == NULL)
    return NULL;

  glGetIntegerv (GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv (GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv (GL_UNPACK_ROW_LENGTH, &rowlength);

  glGetIntegerv (GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv (GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv (GL_UNPACK_ALIGNMENT, &alignment);

  glPixelStorei (GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei (GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei (GL_UNPACK_ROW_LENGTH, 0);

  glPixelStorei (GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei (GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

  glReadPixels (0, 0, (GLsizei)width, (GLsizei)height, format, GL_UNSIGNED_BYTE, (GLvoid*) buffer);

  glPixelStorei (GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei (GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei (GL_UNPACK_ROW_LENGTH, rowlength);
  
  glPixelStorei (GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei (GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei (GL_UNPACK_ALIGNMENT, alignment);
  
  return buffer;
}

#endif
