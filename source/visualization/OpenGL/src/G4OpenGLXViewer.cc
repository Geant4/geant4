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
// $Id: G4OpenGLXViewer.cc 103926 2017-05-03 13:43:27Z gcosmo $
//
// 
// Andrew Walkden  7th February 1997
// G4OpenGLXViewer : Class to provide XWindows specific
//                 functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLXViewer.hh"

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLFontBaseStore.hh"

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4StateManager.hh"
#include "G4VisManager.hh"
#include "G4Text.hh"
#include "G4Threading.hh"

#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <X11/Xmu/StdCmap.h>

#include <assert.h>
#include <sstream>
#include <chrono>
#include <thread>

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
  ((str) != 0 ? (strncpy((char*)malloc((unsigned)strlen(str) + 1), str, (unsigned)strlen(str) + 1)) : (char*)0)

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
#ifdef G4MULTITHREADED
  if (G4Threading::IsMasterThread()) {
    glXMakeCurrent (dpy, win, cxMaster);
  } else {
    glXMakeCurrent (dpy, win, cxVisSubThread);
  }
#else
  glXMakeCurrent (dpy, win, cxMaster);
#endif
  G4OpenGLViewer::SetView ();
}

void G4OpenGLXViewer::ShowView () {
#ifdef G4MULTITHREADED
//  G4int thread_id = G4Threading::G4GetThreadId();
//  G4cout << "G4OpenGLXViewer::ShowView: thread " << thread_id << G4endl;
#endif
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
	  G4cout << Pick(event.xbutton.x, event.xbutton.y) << G4endl;
	}
	else if (event.type == ButtonPress && event.xbutton.button == 2) break;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
}

#ifdef G4MULTITHREADED

void G4OpenGLXViewer::SwitchToVisSubThread()
{
//  G4cout << "G4OpenGLXViewer::SwitchToVisSubThread" << G4endl;
  cxVisSubThread = glXCreateContext (dpy, vi, cxMaster, true);
  glXMakeCurrent (dpy, win, cxVisSubThread);
}

void G4OpenGLXViewer::SwitchToMasterThread()
{
//  G4cout << "G4OpenGLXViewer::SwitchToMasterThread" << G4endl;
  glXMakeCurrent (dpy, win, cxMaster);
  // and destroy sub-thread context
  glXDestroyContext (dpy, cxVisSubThread);
}

#endif

void G4OpenGLXViewer::GetXConnection () {
// get a connection.
  dpy = XOpenDisplay (0);  // Uses DISPLAY environment variable.
  if (!dpy) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXViewer::G4OpenGLXViewer couldn't open display." << G4endl;
    return;
  }

// make sure OpenGL is supported and installed properly.
  if (!glXQueryExtension (dpy, &errorBase, &eventBase)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXViewer::G4OpenGLXViewer X Server has no GLX extension." 
	 << G4endl;
    return;
  }

}

void G4OpenGLXViewer::CreateGLXContext (XVisualInfo* v) {

  vi = v;
// get window's attributes
  if (!XGetWindowAttributes(dpy, XRootWindow (dpy, vi -> screen), &xwa)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXViewer::G4OpenGLXViewer couldn't return window attributes"
	 << G4endl;
    return;
  }
  
// create the master GLX context
  cxMaster = glXCreateContext (dpy, vi, 0, true);
  if (!cxMaster) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXViewer::G4OpenGLXViewer couldn't create context."
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
				      XA_RGB_BEST_MAP,
				      False, 
				      True);

  if (status == 1) {
    cmap = 0;
    status = XGetRGBColormaps (dpy, 
			       XRootWindow (dpy, vi -> screen),
			       &standardCmaps, 
			       &numCmaps, 
			       XA_RGB_BEST_MAP);
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
      if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cerr <<
  "G4OpenGLXViewer::G4OpenGLXViewer failed to allocate a standard colormap."
	       << G4endl;
      return;
    }
    if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "Got standard cmap" << G4endl;
  } else {
    cmap = XCreateColormap (dpy, 
			    XRootWindow(dpy, vi -> screen), 
			    vi -> visual, 
			    AllocNone);
    if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
      G4cout << "Created own cmap" << G4endl;
  }

  if (!cmap) {
    fViewId = -1;  // This flags an error.
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "G4OpenGLXViewer::G4OpenGLXViewer failed to allocate a Colormap."
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
  size_hints = XAllocSizeHints();
    
  ResizeWindow(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());

  G4int x_origin = fVP.GetWindowAbsoluteLocationHintX(DisplayWidth(dpy, vi -> screen));

  // FIXME,  screen size != window size on MAC, but I don't know have to get the menuBar
  // size on MAC. L.Garnier 01/2009
  G4int y_origin = fVP.GetWindowAbsoluteLocationHintY(DisplayHeight(dpy, vi -> screen));

  size_hints->base_width = getWinWidth();
  size_hints->base_height = getWinHeight();
  size_hints->x = x_origin;
  size_hints->y = y_origin;
  if (fVP.IsWindowSizeHintX () && fVP.IsWindowLocationHintX () && fVP.IsWindowLocationHintY ()) {
    size_hints->flags |= PSize | PPosition;
  } else if (fVP.IsWindowSizeHintX () && !(fVP.IsWindowLocationHintX () || fVP.IsWindowLocationHintY ())) {
    size_hints->flags |= PSize;
  } else if ((!fVP.IsWindowSizeHintX ()) && fVP.IsWindowLocationHintX () && fVP.IsWindowLocationHintY ()) {
    size_hints->flags |= PPosition;
  }
  if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "Window name: " << fName << G4endl;
  strncpy (charViewName, fName, 99); charViewName[99] = '\0';
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
                        y_origin, getWinWidth(), getWinHeight(), 0, vi -> depth,
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
  Bool success = glXMakeCurrent (dpy, win, cxMaster);
  if (!success) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXViewer::G4OpenGLXViewer failed to attach a GLX context."
	 << G4endl;
    GLint error = GL_NO_ERROR;
    while ((error = glGetError()) != GL_NO_ERROR) {
      switch (error) {
      case GL_INVALID_ENUM :
	G4cout << "GL Error: GL_INVALID_ENUM" << G4endl;break;
      case GL_INVALID_VALUE :
	G4cout << "GL Error: GL_INVALID_VALUE" << G4endl;break;
      case GL_INVALID_OPERATION :
	G4cout << "GL Error: GL_INVALID_OPERATION" << G4endl;break;
      case GL_OUT_OF_MEMORY :
	G4cout << "GL Error: GL_OUT_OF_MEMORY" << G4endl;break;
      case GL_STACK_UNDERFLOW :
	G4cout << "GL Error: GL_STACK_UNDERFLOW" << G4endl;break;
      case GL_STACK_OVERFLOW :
	G4cout << "GL Error: GL_STACK_OVERFLOW" << G4endl;break;
      default :
	G4cout << "GL Error: " << error << G4endl;break;
      }
    }
    return;
  }
}

void G4OpenGLXViewer::CreateFontLists()
{
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
	"G4OpenGLXViewer::CreateFontLists XLoadQueryFont failed for font\n  "
	     << i->second
	     << G4endl;
      continue;
    }
    G4int font_base = glGenLists(256);
    if (!font_base) {
      G4cerr <<
	"G4OpenGLXViewer::CreateFontLists out of display lists for fonts." 
	     << G4endl;
      continue;
    }
    G4int first = font_info->min_char_or_byte2;
    G4int last  = font_info->max_char_or_byte2;
    glXUseXFont(font_info->fid, first, last-first+1, font_base + first);
    G4int width = font_info->max_bounds.width;
    G4OpenGLFontBaseStore::AddFontBase
      (this, font_base, i->first, i->second, width);
  }
}

void G4OpenGLXViewer::DrawText(const G4Text& g4text)
{
  if (isGl2psWriting()) {

    G4OpenGLViewer::DrawText(g4text);

  } else {

    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text,sizeType);

    const G4OpenGLFontBaseStore::FontInfo& fontInfo =
      G4OpenGLFontBaseStore::GetFontInfo(this,(int)size);
    if (fontInfo.fFontBase < 0) {
      static G4int callCount = 0;
      ++callCount;
      //if (callCount <= 10 || callCount%100 == 0) {
      if (callCount <= 1) {
        G4cout <<
	  "G4OpenGLXViewer::DrawText: No fonts available for \""
	       << fName <<
          "\"\n  Called with "
               << g4text
               << G4endl;
      }
      return;
    }

    const G4Colour& c = fSceneHandler.GetTextColour(g4text);
    glColor4d(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());

    G4Point3D position = g4text.GetPosition();

    G4String textString = g4text.GetText();
    const char* textCString = textString.c_str();
  
    // Set position for raster-style drawers (X, Xm)
    glRasterPos3d(position.x(),position.y(),position.z());

    glPushAttrib(GL_LIST_BIT);

    // Calculate move for centre and right adjustment
    G4double span = textString.size() * fontInfo.fWidth;
    G4double xmove = 0., ymove = 0.;
    switch (g4text.GetLayout()) {
    case G4Text::left: break;
    case G4Text::centre: xmove -= span / 2.; break;
    case G4Text::right: xmove -= span;
    }

    //Add offsets
    xmove += g4text.GetXOffset();
    ymove += g4text.GetYOffset();

    // Do move
    glBitmap(0,0,0,0,xmove,ymove,0);

    // Write characters
    glListBase(fontInfo.fFontBase);
    glCallLists(strlen(textCString),GL_UNSIGNED_BYTE,(GLubyte*)textCString);
    glPopAttrib();
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
  // To satisfy Coverity
  xwa.visual = 0;
  iconName.value = 0;
  xwa.screen = 0;
  windowName.value = 0;

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
    glXDestroyContext (dpy, cxMaster);
    if (win) XDestroyWindow (dpy, win); // ...if already deleted in
    // sub-class G4OpenGLXmViewer.
    XFlush (dpy);
  }
}
	

#endif
