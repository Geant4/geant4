// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXViewer.hh,v 1.11 2001-05-23 14:47:05 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// G4OpenGLXViewer : Class to provide XWindows specific
//                   functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLXVIEWER_HH
#define G4OPENGLXVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "globals.hh"

#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/Xmu/StdCmap.h>

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

class G4OpenGLSceneHandler;

class G4OpenGLXViewer: virtual public G4OpenGLViewer {

public:
  G4OpenGLXViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLXViewer ();
  void SetView ();
  void ShowView ();
  void FinishView ();
  void print();

protected:
  void GetXConnection ();
  void CreateGLXContext (XVisualInfo* vi);
  virtual void CreateMainWindow ();

  char                              print_string[50];
  G4bool                            print_colour,
                                    vectored_ps;

//////////////////////////////Vectored PostScript production functions//////////////////////////////
  void printBuffer(GLint, GLfloat*);
  GLfloat* spewPrimitiveEPS (FILE*, GLfloat*);
  void spewSortedFeedback (FILE*, GLint, GLfloat*);
  void spewWireframeEPS (FILE*, GLint, GLfloat*, const char*);
  void print3DcolorVertex(GLint, GLint*, GLfloat*);
  G4float                           pointSize;


//////////////////////////////Pixmap (screen dump) production functions//////////////////////////////
  GLubyte* grabPixels (int inColor,
		       unsigned int width,
		       unsigned int height);
  int generateEPS (char* filnam,
		   int inColour,
		   unsigned int width,
		   unsigned int height);
  GLXContext create_GL_print_context (XVisualInfo*& pvi, G4bool& db);

  XWindowAttributes                 xwa;
  Display                           *dpy;
  static XVisualInfo                *vi_single_buffer;
  static XVisualInfo                *vi_double_buffer;
  XVisualInfo                       *vi_immediate,
                                    *vi_stored,
                                    *vi;
  Colormap                          cmap;
  XSetWindowAttributes              swa;
  GLXDrawable                       win;
  GLXContext                        cx;
  XEvent                            event;
  G4int                             *attributeList,
                                    errorBase,
                                    eventBase,
                                    major,
                                    minor,
                                    x_origin,
                                    y_origin;
  XSizeHints                        *norm_hints;
  XWMHints                          *wm_hints;
  XClassHint                        *class_hints;
  Pixmap                            icon_pixmap;
  XSizeHints                        *size_hints;
  G4int                             WinSize_x,
                                    WinSize_y;
  Atom                              Xatom;
  XTextProperty                     windowName,
                                    iconName;
  char                              charViewName [100];

private:
  G4OpenGLXViewer (const G4OpenGLXViewer&);
  G4OpenGLXViewer& operator = (const G4OpenGLXViewer&);
};

typedef struct _Feedback3Dcolor {
  GLfloat x;
  GLfloat y;
  GLfloat z;
  GLfloat red;
  GLfloat green;
  GLfloat blue;
  GLfloat alpha;
} Feedback3Dcolor;

#endif

#endif
