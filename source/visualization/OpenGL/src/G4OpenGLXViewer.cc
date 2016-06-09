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
// $Id: G4OpenGLXViewer.cc,v 1.42 2007/05/25 10:47:17 allison Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
// Andrew Walkden  7th February 1997
// G4OpenGLXViewer : Class to provide XWindows specific
//                 functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLXViewer.hh"

#include "G4OpenGLFontBaseStore.hh"

#include <sstream>

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

#include <X11/Xatom.h>
#include <X11/Xutil.h>

#include <assert.h>

int G4OpenGLXViewer::snglBuf_RGBA[12] =
{ GLX_RGBA,
  GLX_RED_SIZE, 1,
  GLX_GREEN_SIZE, 1, 
  GLX_BLUE_SIZE, 1,
  GLX_DEPTH_SIZE, 1,
  GLX_STENCIL_SIZE, 1,
  None };

int G4OpenGLXViewer::dblBuf_RGBA[13] =
{ GLX_RGBA,
  GLX_RED_SIZE, 1,
  GLX_GREEN_SIZE, 1,
  GLX_BLUE_SIZE, 1,
  GLX_DOUBLEBUFFER,
  GLX_DEPTH_SIZE, 1,
  GLX_STENCIL_SIZE, 1,
  None };

#define NewString(str) \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)

#define USE_DEFAULT_COLORMAP 1
#define USE_STANDARD_COLORMAP 0

XVisualInfo*  G4OpenGLXViewer::vi_single_buffer = 0;
XVisualInfo*  G4OpenGLXViewer::vi_double_buffer = 0;

extern "C" {
  static Bool G4OpenGLXViewerWaitForNotify (Display*, XEvent* e, char* arg) {
    return (e->type == MapNotify) && (e->xmap.window == (Window) arg);
  }
}

void G4OpenGLXViewer::SetView () {
  glXMakeCurrent (dpy, win, cx);
  G4OpenGLViewer::SetView ();  
}

void G4OpenGLXViewer::ShowView () {
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propagated before progressing.
  glFlush ();

  if (fVP.IsPicking()) {
    G4cout <<
      "Window activated for picking (left-mouse), exit (middle-mouse)."
	   << G4endl;
    while (true) {
      if (XPending(dpy)) {
	XNextEvent(dpy, &event);
	if (event.type == ButtonPress && event.xbutton.button == 1) {
	  Pick(event.xbutton.x, event.xbutton.y);
	}
	else if (event.type == ButtonPress && event.xbutton.button == 2) break;
      }
    }
  }
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
    cmap = 0;
    status = XGetRGBColormaps (dpy, 
			       XRootWindow (dpy, vi -> screen), 
			       &standardCmaps, 
			       &numCmaps, 
			       XA_RGB_DEFAULT_MAP);
    if (status == 1)
      for (i = 0; i < numCmaps; i++) {
	if (standardCmaps[i].visualid == vi -> visualid) {
	  cmap = standardCmaps[i].colormap;
	  XFree (standardCmaps);
	  break;
	}
      }
    if (!cmap) {
      fViewId = -1;  // This flags an error.
      G4cerr <<
   "G4OpenGLViewer::G4OpenGLViewer failed to allocate a standard colormap."
	     << G4endl;
      return;
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

  // Window size and position...
  unsigned int width, height;
  x_origin = 0;
  y_origin = 0;
  size_hints = XAllocSizeHints();
  const G4String& XGeometryString = fVP.GetXGeometryString();
  int screen_num = DefaultScreen(dpy);
  if (!XGeometryString.empty()) {
    G4int geometryResultMask = XParseGeometry
      ((char*)XGeometryString.c_str(),
       &x_origin, &y_origin, &width, &height);
    if (geometryResultMask & (WidthValue | HeightValue)) {
      if (geometryResultMask & XValue) {
	if (geometryResultMask & XNegative) {
	  x_origin = DisplayWidth(dpy, screen_num) + x_origin - width;
	}
	size_hints->flags |= PPosition;
	size_hints->x = x_origin;
      }
      if (geometryResultMask & YValue) {
	if (geometryResultMask & YNegative) {
	  y_origin = DisplayHeight(dpy, screen_num) + y_origin - height;
	}
	size_hints->flags |= PPosition;
	size_hints->y = y_origin;
      }
    } else {
      G4cout << "ERROR: Geometry string \""
	     << XGeometryString
	     << "\" invalid.  Using \"600x600\"."
	     << G4endl;
      width = 600;
      height = 600;
    }
  }
  size_hints->width = width;
  size_hints->height = height;
  size_hints->flags |= PSize;

  //  G4int                             WinSize_x;
  //  G4int                             WinSize_y;
  WinSize_x = width;
  WinSize_y = height;

  G4cout << "Window name: " << fName << G4endl;
  strncpy (charViewName, fName, 100);
  char *window_name = charViewName;
  char *icon_name = charViewName;
  //char tmpatom[] = "XA_WM_NORMAL_HINTS"; 
  wm_hints = XAllocWMHints();
  class_hints = XAllocClassHint();

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
  XIfEvent (dpy, &event, G4OpenGLXViewerWaitForNotify, (char*) win);

// connect the context to a window
  Bool success = glXMakeCurrent (dpy, win, cx);
  if (!success) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLViewer::G4OpenGLViewer failed to attach a GLX context."
	 << G4endl;
    GLint error = GL_NO_ERROR;
    while ((error = glGetError()) != GL_NO_ERROR) {
      G4cout << "GL Error: " << gluErrorString(error) << G4endl;
    }
    return;
  }

}

void G4OpenGLXViewer::CreateFontLists () {

  std::map<G4double,G4String> fonts;  // G4VMarker screen size and font name.
  fonts[10.] = "-adobe-courier-bold-r-normal--10-100-75-75-m-60-iso8859-1";
  fonts[11.] = "-adobe-courier-bold-r-normal--11-80-100-100-m-60-iso8859-1";
  fonts[12.] = "-adobe-courier-bold-r-normal--12-120-75-75-m-70-iso8859-1";
  fonts[13.] = "fixed";
  fonts[14.] = "-adobe-courier-bold-r-normal--14-100-100-100-m-90-iso8859-1";
  fonts[17.] = "-adobe-courier-bold-r-normal--17-120-100-100-m-100-iso8859-1";
  fonts[18.] = "-adobe-courier-bold-r-normal--18-180-75-75-m-110-iso8859-1";
  fonts[20.] = "-adobe-courier-bold-r-normal--20-140-100-100-m-110-iso8859-1";
  fonts[24.] = "-adobe-courier-bold-r-normal--24-240-75-75-m-150-iso8859-1";
  fonts[25.] = "-adobe-courier-bold-r-normal--25-180-100-100-m-150-iso8859-1";
  fonts[34.] = "-adobe-courier-bold-r-normal--34-240-100-100-m-200-iso8859-1";
  std::map<G4double,G4String>::const_iterator i;
  for (i = fonts.begin(); i != fonts.end(); ++i) {
    XFontStruct* font_info = XLoadQueryFont(dpy, i->second);
    if (!font_info) {
      G4cerr <<
	"G4OpenGLXViewer: XLoadQueryFont failed for font\n  "
	     << i->second
	     << G4endl;
      continue;
    }
    G4int font_base = glGenLists(256);
    if (!font_base) {
      G4cerr << "G4OpenGLXViewer: out of display lists for fonts." 
	     << G4endl;
      continue;
    }
    G4int first = font_info->min_char_or_byte2;
    G4int last  = font_info->max_char_or_byte2;
    glXUseXFont(font_info->fid, first, last-first+1,font_base+first);
    G4OpenGLFontBaseStore::AddFontBase(this,font_base,i->first,i->second);
  }
}

G4OpenGLXViewer::G4OpenGLXViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene),
vi_immediate (0),
vi_stored (0),
vi (0),
cmap (0)
{
  GetXConnection ();
  if (fViewId < 0) return;
  
  // Try for a visual suitable for OpenGLImmediate..
  // first try for a single buffered RGB window
  if (!vi_single_buffer) {
    vi_single_buffer =
      glXChooseVisual (dpy, XDefaultScreen (dpy), snglBuf_RGBA);
  }
  if (!vi_double_buffer) {
    vi_double_buffer =
      glXChooseVisual (dpy, XDefaultScreen (dpy), dblBuf_RGBA);
  }

  if (vi_single_buffer || vi_double_buffer) {
    if (!vi_double_buffer) {
      G4cout <<
	"G4OpenGLXViewer::G4OpenGLXViewer: unable to get a double buffer visual."
	"\n  Working with a single buffer."
	     << G4endl;
    }
  } else {
    if (!vi_single_buffer) {
      G4cout <<
	"G4OpenGLXViewer::G4OpenGLXViewer: unable to get a single buffer visual."
	     << G4endl;
    }
    if (!vi_double_buffer) {
      G4cout <<
	"G4OpenGLXViewer::G4OpenGLXViewer: unable to get a double buffer visual."
	     << G4endl;
    }
  }

  if (vi_single_buffer) {
    vi_immediate = vi_single_buffer;
    attributeList = snglBuf_RGBA;
  }
  
  if (!vi_immediate){
    // next try for a double buffered RGB, but Draw to top buffer
    if (vi_double_buffer) {
      vi_immediate = vi_double_buffer;
      attributeList = dblBuf_RGBA;
    }
  }

  // Now try for a visual suitable for OpenGLStored...
  // Try for a double buffered RGB window
  if (vi_double_buffer) {
    vi_stored = vi_double_buffer;
    attributeList = dblBuf_RGBA;
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
    glXMakeCurrent (dpy, None, NULL);
    glXDestroyContext (dpy, cx);
    if (win) XDestroyWindow (dpy, win); // ...if already deleted in
    // sub-class G4OpenGLXmViewer.
    XFlush (dpy);
  }
}

void G4OpenGLXViewer::print() {
  
  //using namespace std;
  //cout << "print_col_callback requested with file name: " << print_string << G4endl;
  
  if (vectored_ps) {

    G4OpenGLViewer::print();

  } else {

    XVisualInfo* pvi;
    GLXContext pcx = create_GL_print_context(pvi);

    if (!pcx) {
      G4cout << "Unable to create print context." << G4endl;
      return;
    }

    GLXContext tmp_cx;
    tmp_cx = cx;
    cx=pcx;
    
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
		    win,
		    cx);
    
    glViewport (0, 0, WinSize_x, WinSize_y);
    
    ClearView ();
    SetView ();
    DrawView ();
    
    generateEPS (print_string,
		 print_colour,
		 WinSize_x, WinSize_y);
    
    win=tmp_win;
    cx=tmp_cx;
    
    glXMakeCurrent (dpy,
		    win,
		    cx);
    
  }

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

GLXContext G4OpenGLXViewer::create_GL_print_context(XVisualInfo*& pvi) {
  
  pvi = glXChooseVisual (dpy,
			 XDefaultScreen (dpy),
			 snglBuf_RGBA);

  if (!pvi) {
    pvi = glXChooseVisual (dpy,
			   XDefaultScreen (dpy),
			   dblBuf_RGBA);
  }

  return glXCreateContext (dpy,
			   pvi,
			   NULL,
			   False);
}

#endif
